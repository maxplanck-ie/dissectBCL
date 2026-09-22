import json
import logging
import mimetypes
import re
import subprocess as sp
from pathlib import Path
from zipfile import ZIP_STORED, ZipFile

import click
import requests
from rich import print


def _fetch_ro_crate_metadata(request_id, config, parkour_url=None):
    """Fetch RO-Crate metadata from parkour API.

    Args:
        request_id: Request ID extracted from project name
        config: dissectBCL config dict
        parkour_url: Override parkour URL (e.g. for parkour-test)
    """
    url = (parkour_url or config["parkour"]["URL"]).rstrip("/")
    url += "/api/generate_ro_crate/"
    try:
        response = requests.get(
            url,
            params={"requests": request_id, "preview": "true"},
            auth=(config["parkour"]["user"], config["parkour"]["password"]),
            verify=config["parkour"]["cert"],
            timeout=60,
        )
        response.raise_for_status()
        return response.json()["ro_crate"]
    except Exception as e:
        logging.warning(
            f"RO-Crate metadata fetch from {url} failed for request "
            f"{request_id}: {e}. Shipping without ro-crate-metadata.json."
        )
        return None


def _add_fastq_file_entities(ro_crate_metadata, project_dir):
    """Add FASTQ file entities to ro-crate metadata graph."""
    project_dir = Path(project_dir)
    md5sums = {}
    md5_path = project_dir / "md5sums.txt"
    if md5_path.exists():
        for line in md5_path.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) == 2:
                md5sums[parts[0]] = parts[1]

    entities_by_id = {
        entity.get("@id"): entity
        for entity in ro_crate_metadata.get("@graph", [])
        if isinstance(entity, dict)
    }

    for sample_dir in sorted(project_dir.glob("Sample_*")):
        if not sample_dir.is_dir():
            continue
        barcode = sample_dir.name[len("Sample_") :]
        stub_id = f"#fastq-data-{barcode}"
        stub_entity = entities_by_id.get(stub_id)
        if stub_entity is None:
            logging.warning(
                f"RO-Crate: no {stub_id} stub entity found in Parkour metadata; "
                f"skipping fastq file entities for {sample_dir.name}."
            )
            continue

        for fastq_file in sorted(sample_dir.glob("*.fastq.gz")):
            file_id = f"#fastq-file-{barcode}-{fastq_file.name}"
            content_url = f"{project_dir.name}/{sample_dir.name}/{fastq_file.name}"
            encoding_format, _ = mimetypes.guess_type(fastq_file.name)
            file_entity = {
                "@id": file_id,
                "@type": ["File", "MediaObject"],
                "name": fastq_file.name,
                "contentUrl": content_url,
                "encodingFormat": encoding_format or "application/gzip",
            }
            md5 = md5sums.get(fastq_file.name)
            if md5:
                property_id = f"{file_id}-md5"
                file_entity["additionalProperty"] = [{"@id": property_id}]
                ro_crate_metadata["@graph"].append(
                    {
                        "@id": property_id,
                        "@type": "PropertyValue",
                        "name": "md5",
                        "value": md5,
                    }
                )
            ro_crate_metadata["@graph"].append(file_entity)
            entities_by_id[file_id] = file_entity

            has_part = stub_entity.setdefault("hasPart", [])
            file_ref = {"@id": file_id}
            if file_ref not in has_part:
                has_part.append(file_ref)

    return ro_crate_metadata


def _build_ro_crate_archive(project_dir, ro_crate_metadata, fileobj):
    """Build RO-Crate zip and stream to fileobj (e.g. fexsend stdin)."""
    project_dir = Path(project_dir)

    if ro_crate_metadata is not None:
        _add_fastq_file_entities(ro_crate_metadata, project_dir)

    with ZipFile(fileobj, "w", compression=ZIP_STORED) as zip_file:
        # Add all project files
        for file_path in project_dir.rglob("*"):
            if file_path.is_file():
                zip_file.write(
                    file_path, arcname=file_path.relative_to(project_dir.parent)
                )

        # Add ro-crate-metadata.json
        if ro_crate_metadata is not None:
            zip_file.writestr(
                "ro-crate-metadata.json",
                json.dumps(ro_crate_metadata, indent=2),
            )

    # Ensure all data is flushed to the pipe before closing stdin
    fileobj.flush()


def fex(project_path, config, from_address, parkour_url=None):
    """Upload a dissectBCL project to FEX as an RO-Crate archive.

    Args:
        project_path: Path to Project_XXXX_User_PI directory
        config: dissectBCL config dict
        from_address: FEX sender address (from config)
        parkour_url: Override parkour URL (default: use parkour-test for latest fixes)

    The project name must match Project_{request_id}_User_PI format.
    Fetches comprehensive ISA-profile metadata from parkour API and enriches
    it with actual FASTQ file references before streaming to FEX.
    """
    project_dir = Path(project_path)
    if not project_dir.is_dir():
        print(f"[red]Not a directory: {project_path}[/red]")
        return

    project_name = project_dir.name

    # Extract request ID from project name (e.g., Project_3358_Hohl_Manke -> 3358)
    match = re.match(r"^Project_(\d+)_", project_name)
    if not match:
        print(
            f"[red]Cannot extract request ID from project name: {project_name}[/red]\n"
            "[yellow]Expected format: Project_XXXX_User_PI[/yellow]"
        )
        return

    request_id = match.group(1)
    archive_name = f"{project_name}_rocrate.zip"

    # Default to parkour-test for latest ro-crate generation fixes
    if parkour_url is None:
        parkour_url = "https://parkour-test.ie-freiburg.mpg.de"
        print("[dim]Using parkour-test for latest RO-Crate generation code[/dim]")

    print(f"Fetching RO-Crate metadata for request {request_id}...")
    ro_crate_metadata = _fetch_ro_crate_metadata(request_id, config, parkour_url)

    if ro_crate_metadata is None:
        print("[yellow]Proceeding without RO-Crate metadata[/yellow]")

    print(f"Streaming {archive_name} to FEX...")
    fex_proc = sp.Popen(
        ["/home/pipegrp/.local/bin/fexsend", "-s", archive_name, from_address],
        stdin=sp.PIPE,
    )

    try:
        _build_ro_crate_archive(project_dir, ro_crate_metadata, fex_proc.stdin)
    except Exception:
        # Kill fexsend to prevent uploading a truncated/corrupt archive
        fex_proc.kill()
        raise
    finally:
        fex_proc.stdin.close()
        exit_code = fex_proc.wait()

    if exit_code == 0:
        print(f"[green]✓ Uploaded {archive_name} to FEX[/green]")
    else:
        print(f"[red]✗ fexsend exited with code {exit_code}[/red]")
        raise click.Abort()
