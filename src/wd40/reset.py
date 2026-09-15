import shutil
from pathlib import Path

import click
from rich import print

# Illumina demux/QC output that bcl-convert + postmux regenerate on rerun.
ILLUMINA_TARGETS = [
    "Reports",
    "Logs",
    "Project_*",
    "FASTQC_Project_*",
    "Undetermined_*.fastq.gz",
]

# Aviti bases2fastq writes its native output straight into the outLane root
# (not under Reports/), so it needs its own target list.
AVITI_TARGETS = [
    "Samples",
    "info",
    "RunManifest.json",
    "RunParameters.json",
    "IndexAssignment.csv",
    "UnassignedSequences.csv",
    "Metrics.csv",
    "multiqc_report.html",
    "multiqc_data",
    "Logs",
    "Project_*",
    "FASTQC_Project_*",
    "RunStats.json",
    # bases2fastq's own manifest copies, distinct from manifest/RunManifest.csv
    # which dissectBCL reads from and this reset keeps.
    "RunManifest.csv",
]

# Non-hidden and hidden (per-project) done-flags, e.g. bases2fastq.done,
# bclconvert.done, communication.done, .{project}.renamed.done,
# .{project}.postmux.done.
FLAG_GLOBS = ["*.done", ".*.done"]
# fastq.made is the flag a human touches by hand once demux looks good, to
# signal BigRedButton the flowcell is ready to pick up. A reset must clear it
# too, or BRB can start on a lane whose demux output was just wiped.
FLAG_FILES = ["run.failed", "fastq.made"]


def detect_mode(outLane):
    """Aviti keeps its manifest under manifest/RunManifest.csv, Illumina at
    the outLane root as demuxSheet.csv. Returns None if neither is found."""
    outLane = Path(outLane)
    if (outLane / "manifest" / "RunManifest.csv").exists():
        return "aviti"
    if (outLane / "demuxSheet.csv").exists():
        return "illumina"
    return None


def collect_targets(outLane):
    outLane = Path(outLane)
    mode = detect_mode(outLane)
    if mode is None:
        return None, []
    patterns = AVITI_TARGETS if mode == "aviti" else ILLUMINA_TARGETS
    found = []
    for pattern in patterns + FLAG_GLOBS:
        found.extend(outLane.glob(pattern))
    for name in FLAG_FILES:
        candidate = outLane / name
        if candidate.exists():
            found.append(candidate)
    seen = set()
    targets = []
    for path in found:
        if path not in seen:
            seen.add(path)
            targets.append(path)
    return mode, targets


def reset(outLane):
    """Strips an outLane dir back down to its SampleSheet/RunManifest, so it
    can be hand-edited (index mask, mismatches, I5) before a redemux, without
    re-copying from the flowcell's read-only source directory."""
    outLane = Path(outLane)
    mode, targets = collect_targets(outLane)
    if mode is None:
        print(
            f"[red]{outLane} doesn't look like an outLane dir "
            "(no demuxSheet.csv or manifest/RunManifest.csv found). Aborting.[/red]"
        )
        return
    kept = "manifest/RunManifest.csv" if mode == "aviti" else "demuxSheet.csv"
    if not targets:
        print(f"[green]Nothing to reset in {outLane} ({mode} mode).[/green]")
        return
    print(
        f"[bold]About to delete the following {mode} demux output/flags under {outLane}:[/bold]"
    )
    for path in sorted(targets):
        print(f"  {path}")
    print(f"[bold]Keeping:[/bold] {outLane / kept}")
    if not click.confirm("Are you sure?", default=False):
        print("[yellow]Aborted, nothing deleted.[/yellow]")
        return
    for path in targets:
        if path.is_dir() and not path.is_symlink():
            shutil.rmtree(path)
        else:
            path.unlink()
    print(f"[green]Reset complete for {outLane}. Kept {kept}.[/green]")
