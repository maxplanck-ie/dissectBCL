# Contributing to dissectBCL

Thanks for contributing! This project uses automated code-quality checks that
run as [pre-commit] hooks. Please install them before making changes so issues
are caught locally instead of in CI.

## Set up your environment

Install the package with its development dependencies:

```bash
pip install -e .[dev]
```

## Install the git hooks

The hooks are defined in `.pre-commit-config.yaml`. We recommend [`prek`], a
fast drop-in replacement for `pre-commit` (either tool works):

```bash
# Recommended: prek (faster)
prek install

# Or, with pre-commit
pre-commit install
```

Once installed, the hooks run automatically on `git commit` and only against
the files you changed.

## Running the hooks manually

To check the whole repository at any time:

```bash
prek run --all-files      # or: pre-commit run --all-files
```

Most hooks fix issues in place — re-stage the modified files and commit again.

## What the hooks do

| Hook | Purpose |
| --- | --- |
| `trailing-whitespace`, `end-of-file-fixer` | Whitespace hygiene on Python files |
| `check-merge-conflict` | Blocks committing unresolved conflict markers |
| [`codespell`] | Fixes common spelling mistakes in code and docs |
| [`ruff`] (`ruff-check --fix`, `ruff-format`) | Lints and formats Python; config in `pyproject.toml` |
| [`zizmor`] | Security linting for GitHub Actions workflows |

Lint and format rules live under `[tool.ruff]` in `pyproject.toml`. If a
codespell false positive appears (e.g. a variable name), add it to
`ignore-words-list` under `[tool.codespell]` rather than rewording the code.

## Before opening a pull request

- Make sure `prek run --all-files` passes.
- Run the test suite: `pytest`.

[pre-commit]: https://pre-commit.com/
[`prek`]: https://github.com/j178/prek
[`codespell`]: https://github.com/codespell-project/codespell
[`ruff`]: https://docs.astral.sh/ruff/
[`zizmor`]: https://docs.zizmor.sh/
