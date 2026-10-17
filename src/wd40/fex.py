import logging
import re
import shutil
import subprocess as sp
import tempfile
from pathlib import Path

import click
from rich import print

from dissectBCL.ro_crate import (
    STREAMING_MAX_BYTES,
    _fetch_ro_crate_metadata,
    _project_size,
    _sample_barcodes,
    build_ro_crate_archive,
)


# Helpers used here (defined once in dissectBCL.ro_crate, imported above):
# _fetch_ro_crate_metadata, _warn_if_parkour_found_no_records,
# _add_fastq_file_entities, _sample_barcodes, _project_size,
# build_ro_crate_archive, STREAMING_MAX_BYTES
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
            build_ro_crate_archive(
                folders=[project_dir],
                ro_crate_metadata=ro_crate_metadata,
                fileobj=zip_file,
            )
        proc = sp.run([fexsend_path, str(temp_archive), from_address], check=False)
        return proc.returncode
    fex_proc = sp.Popen(
        [fexsend_path, "-s", archive_name, from_address],
        stdin=sp.PIPE,
    )
    try:
        build_ro_crate_archive(
            folders=[project_dir],
            ro_crate_metadata=ro_crate_metadata,
            fileobj=fex_proc.stdin,
        )
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
