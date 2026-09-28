import json
import logging
import mimetypes
import re
import shutil
import subprocess as sp
import tempfile
from pathlib import Path
from zipfile import ZIP_STORED, ZipFile

import click
import requests
from rich import print

# fexsend's streaming (-s) mode is unreliable above this project size: the
# server stores the 56-byte closing MIME boundary in the file and truncates it
# only after the transfer, while fexsend re-checks the stored size ~1 s after
# closing the socket. On large uploads that check sees data+56 and aborts with
# exit 29 (verified on actual 10.6-33 GB projects and a synthetic 10.4 GiB
# repro). Above this threshold we upload the archive as a regular file
# instead, where the server reads an exact byte count and never hits the race.
STREAMING_MAX_BYTES = 2**32


def _fetch_ro_crate_metadata(barcodes, config, parkour_url=None):
    """Fetch RO-Crate metadata from parkour API.

    Args:
        barcodes: Sample/library barcodes (e.g. ["26L004591"]) to look up
        config: dissectBCL config dict
        parkour_url: Override parkour URL (e.g. for parkour-test)

    The parkour ``generate_ro_crate`` endpoint matches the ``barcodes`` query
    parameter against Sample/Library barcodes, and ``requests`` against the
    request *name* (a free-text field, e.g. ``4068_Magat_Oudelaar``). The
    numeric request/pk from a project name like ``Project_4068_...`` is
    neither, so it must be translated to barcodes on the dissectBCL side.
    """
    url = (parkour_url or config["parkour"]["URL"]).rstrip("/")
    url += "/api/generate_ro_crate/"
    try:
        response = requests.get(
            url,
            params={"barcodes": ",".join(barcodes), "preview": "true"},
            auth=(config["parkour"]["user"], config["parkour"]["password"]),
            verify=config["parkour"]["cert"],
            timeout=60,
        )
        response.raise_for_status()
        ro_crate_metadata = response.json()["ro_crate"]
        _warn_if_parkour_found_no_records(ro_crate_metadata, barcodes)
        return ro_crate_metadata
    except Exception as e:
        logging.warning(
            f"RO-Crate metadata fetch from {url} failed for barcodes "
            f"{barcodes}: {e}. Shipping without ro-crate-metadata.json."
        )
        return None


def _warn_if_parkour_found_no_records(ro_crate_metadata, barcodes):
    for entity in ro_crate_metadata.get("@graph", []):
        if (
            entity.get("@id") == "./"
            and entity.get("description")
            == "No matching barcodes or requests were found."
        ):
            logging.warning(
                f"RO-Crate: Parkour returned no records for barcodes "
                f"{', '.join(barcodes)}. fastq file entities will be "
                "synthesized from the project folder instead."
            )
            return


def _add_ref_unique(owner, ref):
    has_part = owner.setdefault("hasPart", [])
    if ref not in has_part:
        has_part.append(ref)


def _synthesize_fastq_data_stub(ro_crate_metadata, entities_by_id, barcode):
    """Create a #fastq-data-<barcode> stub Parkour has not emitted."""
    stub_id = f"#fastq-data-{barcode}"
    stub_entity = {
        "@id": stub_id,
        "@type": "Dataset",
        "name": f"Raw sequencing data for {barcode}",
        "description": (
            "Raw fastq sequencing data, populated at data delivery time by dissectBCL."
        ),
        "identifier": f"urn:parkour:fastq-data:{barcode}",
    }
    ro_crate_metadata["@graph"].append(stub_entity)
    entities_by_id[stub_id] = stub_entity
    root = entities_by_id.get("./")
    if root is not None:
        _add_ref_unique(root, {"@id": stub_id})
    return stub_entity


def _add_fastq_file_entities(ro_crate_metadata, project_dir):
    """Add FASTQ file entities to ro-crate metadata graph.

    Links each Sample_* fastq under its #fastq-data-<barcode> stub. If Parkour
    did not emit the stub (older Parkour versions do not), one is synthesized
    so the delivered RO-Crate still describes every fastq file.
    """
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

    sample_dirs = sorted(
        sample_dir for sample_dir in project_dir.glob("Sample_*") if sample_dir.is_dir()
    )
    missing_stubs = [
        sample_dir
        for sample_dir in sample_dirs
        if f"#fastq-data-{sample_dir.name[len('Sample_') :]}" not in entities_by_id
    ]
    if missing_stubs:
        logging.warning(
            f"RO-Crate: Parkour did not provide {len(missing_stubs)} "
            "#fastq-data-* stub entities; synthesizing them so fastq files "
            "are still described."
        )

    for sample_dir in sample_dirs:
        barcode = sample_dir.name[len("Sample_") :]
        stub_id = f"#fastq-data-{barcode}"
        stub_entity = entities_by_id.get(stub_id)
        if stub_entity is None:
            stub_entity = _synthesize_fastq_data_stub(
                ro_crate_metadata, entities_by_id, barcode
            )

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

            _add_ref_unique(stub_entity, {"@id": file_id})

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


def _sample_barcodes(project_dir):
    """Return sample barcodes from Sample_* subdirectories, in run order."""
    project_dir = Path(project_dir)
    return sorted(
        sample_dir.name[len("Sample_") :]
        for sample_dir in project_dir.glob("Sample_*")
        if sample_dir.is_dir()
    )


def _project_size(project_dir):
    """Total size in bytes of every regular file under project_dir."""
    return sum(p.stat().st_size for p in Path(project_dir).rglob("*") if p.is_file())


def _fex_archive_exists(archive_name, from_address):
    """Return whether archive_name is already present on the FEX server."""
    fex_list = sp.check_output(["fexsend", "-l", from_address]).decode("utf-8")
    return archive_name in fex_list.replace("\n", " ").split(" ")


def _fex_delete(archive_name, from_address):
    """Delete archive_name from the FEX server, if present."""
    if _fex_archive_exists(archive_name, from_address):
        sp.run(["fexsend", "-d", archive_name, from_address], check=False)


def _fex_upload(
    fexsend_path,
    archive_name,
    from_address,
    project_dir,
    ro_crate_metadata,
    temp_archive=None,
):
    """Run one upload attempt; returns the fexsend exit code.

    Streaming mode builds the zip straight into fexsend's stdin (no temp
    copy). For large projects temp_archive is a path: the zip is built there
    first and uploaded as a regular file, which avoids fexsend's streaming
    byte-mismatch race on multi-GiB payloads.
    """
    if temp_archive is not None:
        with open(temp_archive, "wb") as zip_file:
            _build_ro_crate_archive(project_dir, ro_crate_metadata, zip_file)
        proc = sp.run([fexsend_path, str(temp_archive), from_address], check=False)
        return proc.returncode

    fex_proc = sp.Popen(
        [fexsend_path, "-s", archive_name, from_address],
        stdin=sp.PIPE,
    )
    try:
        _build_ro_crate_archive(project_dir, ro_crate_metadata, fex_proc.stdin)
    except Exception:
        # Kill fexsend to prevent uploading a truncated/corrupt archive
        fex_proc.kill()
        _fex_delete(archive_name, from_address)
        raise
    finally:
        fex_proc.stdin.close()
    return fex_proc.wait()


def fex(project_path, config, from_address, parkour_url=None):
    """Upload a dissectBCL project to FEX as an RO-Crate archive.

    Args:
        project_path: Path to Project_XXXX_User_PI directory
        config: dissectBCL config dict
        from_address: FEX sender address (from config)
        parkour_url: Override parkour URL (default: use URL from config)

    The project name must match Project_{request_id}_User_PI format.
    Fetches comprehensive ISA-profile metadata from parkour API and enriches
    it with actual FASTQ file references before streaming to FEX.
    """
    project_dir = Path(project_path)
    if not project_dir.is_dir():
        print(f"[red]Not a directory: {project_path}[/red]")
        return

    project_name = project_dir.name

    # Validate project name matches Project_{request_id}_User_PI format,
    # e.g. Project_3358_Hohl_Manke
    match = re.match(r"^Project_(\d+)_", project_name)
    if not match:
        print(
            f"[red]Cannot extract request ID from project name: {project_name}[/red]\n"
            "[yellow]Expected format: Project_XXXX_User_PI[/yellow]"
        )
        return

    archive_name = f"{project_name}_rocrate.zip"

    # Use parkour URL from config, allow override via --parkour-url
    if parkour_url is None:
        parkour_url = config["parkour"]["URL"]

    barcodes = _sample_barcodes(project_dir)
    if barcodes:
        print(f"Fetching RO-Crate metadata for {len(barcodes)} sample barcodes...")
        ro_crate_metadata = _fetch_ro_crate_metadata(barcodes, config, parkour_url)
    else:
        logging.warning(
            f"RO-Crate: no Sample_* folders found under {project_name}; "
            "shipping without ro-crate-metadata.json."
        )
        ro_crate_metadata = None

    if ro_crate_metadata is None:
        print("[yellow]Proceeding without RO-Crate metadata[/yellow]")

    # Find fexsend: check PATH first, fall back to known location
    fexsend_path = shutil.which("fexsend") or "/home/pipegrp/.local/bin/fexsend"

    # Never stack up stale archives: drop any previous upload of the same name
    # (a fresh archive is built on every run), so only the correct file remains.
    _fex_delete(archive_name, from_address)

    # fexsend's streaming mode mis-verifies large uploads (see
    # STREAMING_MAX_BYTES); for those, build the zip in a temp dir on the same
    # filesystem and upload it as a regular file. The temp archive is removed
    # afterwards, so the delivered copy lives only on the FEX server. The temp
    # dir sits next to the project (same mount) with the archive_name as the
    # file's basename, so fexsend stores it under the expected name.
    project_size = _project_size(project_dir)
    tmp_dir = None
    temp_archive = None
    if project_size >= STREAMING_MAX_BYTES:
        tmp_dir = Path(tempfile.mkdtemp(prefix=".fextmp_", dir=project_dir.parent))
        temp_archive = tmp_dir / archive_name
        print(
            f"Archive is {project_size / 2**30:.1f} GiB; building it to a temp file "
            "and uploading as a regular file (avoids the fexsend streaming race)."
        )
    else:
        print(f"Streaming {archive_name} to FEX...")

    try:
        last_exit_code = None
        for attempt in (1, 2):
            try:
                last_exit_code = _fex_upload(
                    fexsend_path,
                    archive_name,
                    from_address,
                    project_dir,
                    ro_crate_metadata,
                    temp_archive=temp_archive,
                )
            except Exception:
                # A failed build uploads nothing; _fex_upload already cleaned up
                # any streaming fexsend. Re-raise so the caller sees the error.
                raise

            if last_exit_code == 0:
                print(f"[green]✓ Uploaded {archive_name} to FEX[/green]")
                return

            # fexsend detected a mismatch between streamed bytes and what the
            # server stored (usually a shelled 56-byte multipart boundary). The
            # server copy is corrupt, so remove it before retrying to avoid
            # replacing a stale archive with another broken one.
            print(
                f"[red]✗ fexsend exited with code {last_exit_code} "
                f"(attempt {attempt}/2); deleting the corrupt upload and retrying.[/red]"
            )
            _fex_delete(archive_name, from_address)

        print(f"[red]✗ fexsend exited with code {last_exit_code} after retrying.[/red]")
        raise click.Abort()
    finally:
        if tmp_dir is not None:
            shutil.rmtree(tmp_dir, ignore_errors=True)
