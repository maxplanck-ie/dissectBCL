import glob
import json
import logging
import os
import subprocess as sp

import requests
from rich import print


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


def hashFile(fil):
    """
    BLAKE3 of a single file via the `b3sum` binary (conda-forge,
    see env.yml). Multithreaded internally and considerably faster
    than hashlib.md5 on large BAM/fastq files.
    """
    out = sp.run(["b3sum", "--no-names", fil], stdout=sp.PIPE, check=True)
    return out.stdout.decode().strip()


def dirhash(path):
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
            entries.append(f"{rel}:{hashFile(fil)}")
    entries.sort()
    manifest = "\n".join(entries).encode()
    out = sp.run(
        ["b3sum", "--no-names"],
        input=manifest,
        stdout=sp.PIPE,
        check=True,
    )
    return out.stdout.decode().strip()


def pushChecksum(parkourURL, parkourAuth, parkourCert, reqID, path, checksum):
    d = {"data": path, "md5": checksum}
    r = requests.post(
        parkourURL + "/api/requests/" + reqID + "/put_filepaths/",
        auth=parkourAuth,
        data=d,
        verify=parkourCert,
    )
    return r


def drain(configpath, parkourURL, parkourAuth, parkourCert):
    """Process every queued job once. Safe to call repeatedly (e.g. cron)."""
    d = queueDir(configpath)
    for job in sorted(glob.glob(os.path.join(d, "*.json"))):
        with open(job) as fh:
            payload = json.load(fh)
        reqID, path = payload["reqID"], payload["path"]
        if not os.path.exists(path):
            print(f"[red]{path} no longer exists, dropping job.[/red]")
            os.remove(job)
            continue
        print(f"Hashing {path} for reqID {reqID}...")
        checksum = dirhash(path)
        r = pushChecksum(parkourURL, parkourAuth, parkourCert, reqID, path, checksum)
        if r.status_code == 200:
            print(f"[green]Pushed checksum for {reqID}:[/green] {checksum}")
            os.remove(job)
        else:
            print(
                f"[red]Failed pushing checksum for {reqID} "
                f"({r.status_code}), will retry next run.[/red]"
            )
