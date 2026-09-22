# AGENTS.md

Guidance for agents (and humans) working in this repository.

## What this is

dissectBCL demultiplexes Illumina (novaseq/miseq/nextseq) and Aviti sequencing
runs. It is a Python package (`src/dissectBCL`, `src/wd40`, `src/tools`), with
CLI entry points `dissect`, `wd40`, `email`, and `contam` (see `pyproject.toml`).
It integrates with Parkour2, the LIMS, through a REST API: it pulls sample
metadata before demux, and pushes run QC statistics after.

## Environment

See `README.md` for the production install path (`install_dissect.sh`,
which creates a version-named conda environment per release). For the dev
loop, in an existing checkout: create the conda environment from `env.yml`,
then run `pip install .` (or `.[dev]`, for tests and lint). Versioning is
derived from git tags, through `setuptools_scm`, so `git` must be present
in the environment, and the repo must be a real git checkout, not a
tarball with no `.git` folder.

## Code quality

- Pre-commit hooks are configured in `.pre-commit-config.yaml`
  (trailing-whitespace, end-of-file-fixer, check-merge-conflict, codespell,
  ruff check/format, zizmor). Install them with `prek install`, or with
  `pre-commit install`, before you make changes — see `CONTRIBUTING.md` for
  full details.
- Run `prek run --all-files` (or `pre-commit run --all-files`) and `pytest`,
  before you open a PR.

## Commit messages

This repo uses [Conventional Commits](https://www.conventionalcommits.org/)
(`feat:`, `fix:`, `chore:`, `docs:`, and more). `release-please` generates
releases and `CHANGELOG.md` entries (see `release-please-config.json`) by
reading commit messages to decide the version bump. A commit message that
does not follow this format will not be picked up correctly.

## Tests

Tests live in `tests/` (pytest). Fixtures and sample data for demux tests are
under `tests/test_demux/`.

## wd40 fex: Uploading RO-Crate archives to FEX

The `wd40 fex` subcommand uploads a dissectBCL project to FEX as an RO-Crate
archive, following the same pattern as `dissectBCL.misc.fexUpload` but exposed
as a standalone command:

```bash
wd40 fex Project_3358_Hohl_Manke
```

This:
1. Extracts the request ID from the project name (`Project_3358_Hohl_Manke` → `3358`)
2. Fetches comprehensive ISA-profile RO-Crate metadata from parkour API
   (`/api/generate_ro_crate/?requests=3358&preview=true`)
   - Uses URL from config by default
   - Can override with `--parkour-url` flag
3. Adds FASTQ file entities to the metadata:
   - Reads `md5sums.txt` for checksums
   - Creates `#fastq-file-{barcode}-{filename}` entities
   - Links them to `#fastq-data-{barcode}` stubs from parkour metadata
4. Builds a zip archive: all project files + `ro-crate-metadata.json`
5. Streams the zip directly to `fexsend` without writing to disk

**Parkour URL:** Uses `parkour.URL` from config by default (production). Override
with `--parkour-url` to test with parkour-test or parkour-dev:

```bash
wd40 fex --parkour-url https://parkour-test.ie-freiburg.mpg.de Project_3358_Hohl_Manke
```

**Requirements:**
- Project name must match `Project_{request_id}_User_PI` format
- Project must exist in parkour database under that request ID
- `fexsend` must be in PATH (typically `~/.local/bin/fexsend` on rapidus)
- Config file must have parkour credentials (`~/configs/dissectBCL_*.ini`)
