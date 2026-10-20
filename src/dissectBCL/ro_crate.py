"""Shared RO-Crate metadata fetch, augmentation, and archive building for dissectBCL.

Consolidates logic previously duplicated in misc.py and wd40/fex.py.
"""

import json
import logging
import mimetypes
from pathlib import Path
from zipfile import ZIP_STORED, ZipFile

import requests

# fexsend's streaming (-s) mode is unreliable above this project size: the
# server stores the 56-byte closing MIME boundary in the file and truncates it
# only after the transfer, while fexsend re-checks the stored size ~1 s after
# closing the socket. On large uploads that check sees data+56 and aborts with
# exit 29 (verified on actual 10.6-33 GB projects and a synthetic 10.4 GiB
# repro). Above this threshold we upload the archive as a regular file
# instead, where the server reads an exact byte count and never hits the race.
STREAMING_MAX_BYTES = 2**32


def _add_ref_unique(owner, ref):
    """Add a reference to owner.hasPart if not already present."""
    has_part = owner.setdefault("hasPart", [])
    if ref not in has_part:
        has_part.append(ref)


def _fetch_ro_crate_metadata(barcodes, config, parkour_url=None):
    """Fetch RO-Crate metadata from Parkour API.

    Args:
        barcodes: Sample/library barcodes (e.g. ["26L004591"]) to look up
        config: dissectBCL config dict (with "parkour" section)
        parkour_url: Override parkour URL (e.g. for parkour-test)

    The Parkour ``generate_ro_crate`` endpoint matches the ``barcodes`` query
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
            f"RO-Crate metadata fetch from Parkour failed for barcodes "
            f"{barcodes}: {e}. Shipping without ro-crate-metadata.json."
        )
        return None


def _warn_if_parkour_found_no_records(ro_crate_metadata, barcodes):
    """Log a warning if Parkour returned an empty crate for the given barcodes."""
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
            break


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


def _project_size(project_dir):
    """Total size in bytes of every regular file under project_dir."""
    return sum(p.stat().st_size for p in Path(project_dir).rglob("*") if p.is_file())


def _sample_barcodes(project_dir):
    """Return sample barcodes from Sample_* subdirectories, in run order."""
    return sorted(
        sample_dir.name[len("Sample_") :]
        for sample_dir in Path(project_dir).glob("Sample_*")
        if sample_dir.is_dir()
    )


def build_ro_crate_archive(
    folders,
    ro_crate_metadata,
    fileobj=None,
    archive_name=None,
):
    """Build an RO-Crate zip archive.

    Args:
        folders: list of directories to include. Each is zipped with paths
            relative to its own parent, so each folder contributes its own
            basename to the archive (e.g. Project_3358_Hohl_Manke/Sample_x/...).
            folders[0] is the project dir containing Sample_* folders and is
            the one the FASTQ entities are augmented onto.
        ro_crate_metadata: the crate dict; FASTQ file entities are added in
            place. May be None, in which case no metadata file is written.
        fileobj: if given, stream the zip there (e.g. a subprocess's stdin
            pipe), avoiding a second on-disk copy of the fastq data.
        archive_name: name for the on-disk archive when fileobj is None.

    Returns the written Path when fileobj is None, else None.
    """
    folders = [Path(f) for f in folders]

    if ro_crate_metadata is not None:
        _add_fastq_file_entities(ro_crate_metadata, folders[0])

    if fileobj is not None:
        target = fileobj
    else:
        if archive_name is None:
            archive_name = f"{folders[0].name}_ro_crate.zip"
        target = folders[0].parent / archive_name

    with ZipFile(target, "w", compression=ZIP_STORED) as zip_file:
        for folder in folders:
            for file_path in folder.rglob("*"):
                if file_path.is_file():
                    zip_file.write(
                        file_path, arcname=file_path.relative_to(folder.parent)
                    )
        if ro_crate_metadata is not None:
            zip_file.writestr(
                "ro-crate-metadata.json",
                json.dumps(ro_crate_metadata, indent=2),
            )

    return None if fileobj is not None else target
