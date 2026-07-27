import glob
import json
import logging
import os
import subprocess as sp
from pathlib import Path

import requests
from rich import print


class ChecksumInterrupted(Exception):
    """Raised when a dissect instance starts demuxing mid-hash."""


def queueDir(configpath):
    """
    Shared spool directory for pending checksum jobs, one .json file
    per released project. Queued by `wd40 rel`, drained by
    `wd40 checksum` (meant to run off a systemd timer / cron, since
    hashing a whole Analysis_* dir is too slow to do inline in `rel`).
    """
    d = os.path.join(
        os.path.expanduser("~/.wd40_checksum_queue"),
        os.path.basename(os.path.normpath(configpath)),
    )
    os.makedirs(d, exist_ok=True)
    return d


def enqueue(configpath, reqID, path):
    """Called from wd40.release.rel() right after put_filepaths succeeds."""
    d = queueDir(configpath)
    job = os.path.join(d, f"{reqID}__{os.path.basename(path)}.json")
    with open(job, "w") as fh:
        json.dump({"reqID": reqID, "path": path}, fh)
    logging.info(f"checksum - queued {path} for reqID {reqID}")


def isBusy(busyDir):
    """
    True if any dissect instance (aviti, illumina, ...) is currently
    demultiplexing a flowcell - see dissectBCL.misc.busyDir/setBusy/
    clearBusy, which is the single directory all instances and this
    tool agree on regardless of each instance's own ini.
    """
    return any(Path(busyDir).glob("dissect.*.busy"))


def b3sumCmd(numThreads=0):
    """
    Base `b3sum` invocation. b3sum uses all available cores by
    default (rayon's global thread pool), so numThreads=0 (the
    default) just omits the flag. Pass a positive int to cap it, e.g.
    on a shared/constrained host.
    """
    cmd = ["b3sum", "--no-names"]
    if numThreads:
        cmd += ["--num-threads", str(numThreads)]
    return cmd


def runB3sum(cmd, input_bytes=None, busyDir=None, pollInterval=2):
    """
    Run a b3sum command, polling `busyDir` every `pollInterval`
    seconds while it's running. If a dissect instance shows up as
    busy mid-hash, the subprocess is killed and ChecksumInterrupted is
    raised - the whole point of doing this async in the first place is
    to never compete with demux for cores, even on a single huge file.
    """
    proc = sp.Popen(
        cmd,
        stdin=sp.PIPE if input_bytes is not None else sp.DEVNULL,
        stdout=sp.PIPE,
    )
    pending_input = input_bytes
    while True:
        if busyDir and isBusy(busyDir):
            proc.kill()
            proc.wait()
            raise ChecksumInterrupted()
        try:
            out, _ = proc.communicate(input=pending_input, timeout=pollInterval)
            break
        except sp.TimeoutExpired:
            pending_input = None
    if proc.returncode != 0:
        raise sp.CalledProcessError(proc.returncode, cmd)
    return out.decode().strip()


def hashFile(fil, numThreads=0, busyDir=None):
    """
    BLAKE3 of a single file via the `b3sum` binary (conda-forge,
    see env.yml). Multithreaded internally and considerably faster
    than hashlib.md5 on large BAM/fastq files.
    """
    return runB3sum(b3sumCmd(numThreads) + [fil], busyDir=busyDir)


def dirhash(path, numThreads=0, busyDir=None):
    """
    Aggregate hash for a directory: hash every regular file, then
    hash the sorted "relpath:filehash" listing itself. Same shape as
    tools like `dirhash`/`checksumdir`, just backed by b3sum instead
    of md5 for speed.
    """
    entries = []
    for root, _, files in os.walk(path):
        for f in files:
            fil = os.path.join(root, f)
            if os.path.islink(fil):
                continue
            rel = os.path.relpath(fil, path)
            entries.append(f"{rel}:{hashFile(fil, numThreads, busyDir)}")
    entries.sort()
    manifest = "\n".join(entries).encode()
    return runB3sum(b3sumCmd(numThreads), input_bytes=manifest, busyDir=busyDir)


def pushChecksum(parkourURL, parkourAuth, parkourCert, reqID, path, checksum):
    d = {"data": path, "md5": checksum}
    r = requests.post(
        parkourURL + "/api/requests/" + reqID + "/put_filepaths/",
        auth=parkourAuth,
        data=d,
        verify=parkourCert,
    )
    return r


def drain(configpath, parkourURL, parkourAuth, parkourCert, numThreads=0, busyDir=None):
    """Process every queued job once. Safe to call repeatedly (e.g. cron)."""
    d = queueDir(configpath)
    for job in sorted(glob.glob(os.path.join(d, "*.json"))):
        if busyDir and isBusy(busyDir):
            print(
                "[yellow]dissect is busy demultiplexing a flowcell, "
                "leaving remaining checksum jobs queued for next run.[/yellow]"
            )
            break
        with open(job) as fh:
            payload = json.load(fh)
        reqID, path = payload["reqID"], payload["path"]
        if not os.path.exists(path):
            print(f"[red]{path} no longer exists, dropping job.[/red]")
            os.remove(job)
            continue
        print(f"Hashing {path} for reqID {reqID}...")
        try:
            checksum = dirhash(path, numThreads, busyDir)
        except ChecksumInterrupted:
            print(
                f"[yellow]Interrupted hashing {path} - dissect started "
                "demuxing. Job stays queued for next run.[/yellow]"
            )
            break
        r = pushChecksum(parkourURL, parkourAuth, parkourCert, reqID, path, checksum)
        if r.status_code == 200:
            print(f"[green]Pushed checksum for {reqID}:[/green] {checksum}")
            os.remove(job)
        else:
            print(
                f"[red]Failed pushing checksum for {reqID} "
                f"({r.status_code}), will retry next run.[/red]"
            )
