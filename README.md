[![Documentation Status](https://readthedocs.org/projects/dissectbcl/badge/?version=latest)](https://dissectbcl.readthedocs.io/en/latest/?badge=latest)
[![Lint](https://github.com/maxplanck-ie/dissectBCL/actions/workflows/lint.yml/badge.svg)](https://github.com/maxplanck-ie/dissectBCL/actions/workflows/lint.yml)
![Pytest](https://github.com/maxplanck-ie/dissectBCL/actions/workflows/pytest.yml/badge.svg)

# dissectBCL

Demultiplexing pipeline for illumina data (novaseq/miseq/nextseq). Continuation of Devon Ryan's [TWTWTWTW](https://github.com/maxplanck-ie/TheWhoTheWhatTheHuh).

## Installation.

Clone this repository and run the install script. It creates a version-named
conda env (e.g. `dissect_v1.0.3`) from the latest git tag and installs
dissectBCL into it, so the reported version always matches what's running
and older releases stay available side by side.

 > git clone git@github.com:maxplanck-ie/dissectBCL.git  
 > cd dissectBCL  
 > ./install_dissect.sh  
 > conda activate dissect_v1.0.3  

Run `./install_dissect.sh` again after pulling a new release (e.g. once
release-please tags and merges a new version) to create the next versioned
env; it prunes old envs automatically, keeping the 5 most recent. See
`./install_dissect.sh -h` for options (specific tag, env prefix, how many
to keep, force-recreate).

> [!NOTE]
> If you get `LookupError: setuptools-scm was unable to detect version`, ensure git is available in the conda environment: `conda install git`.

## Running.

Fill in the dissectBCL.ini file appropriately. By default the config file is expected to be in ~/configs/dissectBCL_prod.ini.

 > dissect

or 

 > dissect -c /path/to/config.ini

or

 > dissect -f /path/to/flowcell.ini

## wd40 checksum timer.

`wd40 rel` releases a flowcell to the periphery and pushes filepaths to
Parkour2 immediately, but it only *queues* a checksum job rather than
hashing the (potentially huge) released directories inline. A separate
`wd40 checksum` command drains that queue: it hashes each queued
directory with [`b3sum`](https://github.com/BLAKE3-team/BLAKE3)
(BLAKE3, installed via `env.yml`/conda-forge) and pushes the result to
Parkour2 as the `md5` field alongside the existing filepath entry.

Run it periodically via a systemd **user** timer (unit files in
[`systemd/`](systemd/)):

 > mkdir -p ~/.config/systemd/user  
 > cp systemd/wd40-checksum.service systemd/wd40-checksum.timer ~/.config/systemd/user/  

Edit `~/.config/systemd/user/wd40-checksum.service` and point
`ExecStart` at the `wd40` binary inside the conda env you installed
dissectBCL into (see Installation above, e.g.
`~/miniconda3/envs/dissect_v1.0.4/bin/wd40 checksum`).

Then enable and start the timer:

 > systemctl --user daemon-reload  
 > systemctl --user enable --now wd40-checksum.timer  

Check it's scheduled, or run a job manually:

 > systemctl --user list-timers wd40-checksum.timer  
 > systemctl --user start wd40-checksum.service   # run once now  
 > journalctl --user -u wd40-checksum.service      # logs  

If your user session doesn't stay alive after logout (no active
login), enable lingering once so the timer keeps firing:

 > loginctl enable-linger $USER  

### Configuring cores.

`b3sum` uses all available cores by default. To cap it (e.g. a host
shared with other jobs), set `checksumThreads` under `[wd40]` in your
`dissectBCL.ini`:

```ini
[wd40]
checksumThreads=0   # 0 = all available cores (default)
```

Or override it per-run without touching the config:

 > wd40 checksum --threads 8

### Yielding to `dissect`.

If you run multiple `dissect` instances on one host (e.g. one for
Aviti, one for Illumina), they can demux at the same time and `wd40
checksum` needs to stay out of the way of *either* of them.

Every `dissect` instance touches `<busyDir>/dissect.<pid>.busy` for
as long as it's processing a flowcell (`prepConvert` through
`organiseLogs`, cleared in a `finally` so a crash can't leave it
stuck), where `busyDir` is `[wd40] busyDir` from the ini -
`~/.dissectBCL/busy` by default. It's deliberately *not* under each
instance's own `tempDir`: all `dissect`/`wd40` ini files on a host
must point `busyDir` at the same shared location so one `wd40
checksum` can see every instance's flag. PID-suffixed filenames mean
concurrent instances never clobber each other's flag.

`wd40 checksum` only starts a job when that directory is completely
empty of `*.busy` files (i.e. the server is fully idle w.r.t.
demultiplexing, not just "the instance that queued this job"), and
polls it every couple of seconds *while* hashing - including mid-file
on a single large hash - killing `b3sum` and leaving the job queued
for the next timer tick the moment any instance starts demuxing. This
means hashing never competes with demux for cores, regardless of
which sequencer triggered it. On a setup where `busyDir` isn't shared
correctly, the systemd unit's `Nice=`/`IOSchedulingClass=` settings
are the fallback deprioritization.

> [!IMPORTANT]
> This only works if both `dissectBCL_aviti.ini` and
> `dissectBCL_illumina.ini` are updated to point at the same `busyDir`
> value - that's a manual config step on your end when you deploy
> this, not something the code can enforce.

## Docs.

Documentation is available [here](https://dissectbcl.readthedocs.io/en/latest/).
