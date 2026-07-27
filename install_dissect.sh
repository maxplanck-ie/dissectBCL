#!/usr/bin/env bash
#
# install_dissect.sh - create a version-named conda env for dissectBCL and
# install a tagged release into it non-editably. Not part of the Python
# package; a local ops convenience script run after pulling a new release.
#
# Usage: ./install_dissect.sh [-t <tag>] [-p <prefix>] [-k <n>] [-f] [-s] [-h]

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: install_dissect.sh [-t <tag>] [-p <prefix>] [-k <n>] [-f] [-s] [-h]

Create a version-named conda env (e.g. dissect_v1.0.4) and install
dissectBCL into it from a git tag, then prune old versioned envs.

Options:
  -t <tag>     Install this tag instead of the latest reachable tag
               (e.g. v1.0.3). Default: latest tag from `git describe`.
  -p <prefix>  Env name prefix. Default: dissect_v
  -k <n>       Number of versioned envs to keep (newest first). Default: 5
  -f           Force: remove the target env first if it already exists.
  -s           Also (re)install the wd40-checksum systemd --user timer,
               pointed at the wd40 binary of the env just installed.
  -h           Show this help and exit.
EOF
}

prefix="dissect_v"
keep=5
tag=""
force=0
install_systemd=0

while getopts ":t:p:k:fsh" opt; do
  case "$opt" in
    t) tag="$OPTARG" ;;
    p) prefix="$OPTARG" ;;
    k) keep="$OPTARG" ;;
    f) force=1 ;;
    s) install_systemd=1 ;;
    h) usage; exit 0 ;;
    \?) echo "Unknown option: -$OPTARG" >&2; usage; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument" >&2; usage; exit 1 ;;
  esac
done

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

command -v conda >/dev/null 2>&1 || { echo "conda not found on PATH" >&2; exit 1; }
command -v git   >/dev/null 2>&1 || { echo "git not found on PATH" >&2; exit 1; }
[ -d .git ] || { echo "$script_dir is not a git repository" >&2; exit 1; }

git fetch --tags --quiet

if [ -z "$tag" ]; then
  tag="$(git describe --tags --abbrev=0)"
fi

if ! git rev-parse "refs/tags/$tag" >/dev/null 2>&1; then
  echo "Tag '$tag' not found (did you git fetch --tags?)" >&2
  exit 1
fi

version="${tag#v}"
envname="${prefix}${version}"

# Build from a throwaway worktree at the tag rather than checking out this
# repo in place - avoids touching the caller's working tree (which may hold
# this very script, absent from older tags) or requiring it to be clean.
worktree_dir="$(mktemp -d "${TMPDIR:-/tmp}/dissect_install_${version}.XXXXXX")"
cleanup() { git worktree remove --force "$worktree_dir" >/dev/null 2>&1 || rm -rf "$worktree_dir"; }
trap cleanup EXIT

echo "Checking out $tag into $worktree_dir..."
git worktree add --quiet --detach "$worktree_dir" "$tag"

[ -f "$worktree_dir/env.yml" ] || { echo "env.yml not found at $tag" >&2; exit 1; }

env_exists=0
conda env list | awk '{print $1}' | grep -qx "$envname" && env_exists=1

if [ "$force" -eq 1 ]; then
  # Unconditional: env listings can lag behind actual env directories on
  # some shared filesystems, so don't trust env_exists alone to decide
  # whether removal is needed - just make removal a no-op if absent.
  echo "Removing existing env $envname (forced, if present)..."
  conda env remove -n "$envname" --yes >/dev/null 2>&1 || true
elif [ "$env_exists" -eq 1 ]; then
  echo "Env '$envname' already exists. Use -f to remove and recreate it." >&2
  exit 1
fi

echo "Creating env $envname from env.yml..."
conda env create -f "$worktree_dir/env.yml" --name "$envname"

echo "Installing dissectBCL $tag into $envname..."
conda run -n "$envname" pip install "$worktree_dir"

if [ "$install_systemd" -eq 1 ]; then
  if [ ! -f "$worktree_dir/systemd/wd40-checksum.service" ]; then
    echo "WARNING: -s requested but $tag predates systemd/wd40-checksum.service; skipping." >&2
  elif command -v systemctl >/dev/null 2>&1; then
    echo "Installing wd40-checksum systemd --user units..."
    wd40_bin="$(conda run -n "$envname" command -v wd40)"
    systemd_user_dir="${XDG_CONFIG_HOME:-$HOME/.config}/systemd/user"
    mkdir -p "$systemd_user_dir"
    sed "s|@WD40_BIN@|$wd40_bin|" \
      "$worktree_dir/systemd/wd40-checksum.service" \
      > "$systemd_user_dir/wd40-checksum.service"
    cp "$worktree_dir/systemd/wd40-checksum.timer" "$systemd_user_dir/wd40-checksum.timer"
    systemctl --user daemon-reload
    systemctl --user enable --now wd40-checksum.timer
    echo "wd40-checksum.timer enabled (points at $wd40_bin)."
  else
    echo "WARNING: -s requested but systemctl not found; skipping systemd install." >&2
  fi
fi

echo "Pruning old '$prefix*' envs, keeping newest $keep..."
mapfile -t old_envs < <(
  conda env list | awk '{print $1}' \
    | grep -E "^${prefix}[0-9]" \
    | sort -Vr \
    | tail -n +"$((keep + 1))"
)
for old in "${old_envs[@]:-}"; do
  [ -z "$old" ] && continue
  echo "Removing old env: $old"
  conda env remove -n "$old" --yes
done

cat <<EOF

Installed dissectBCL $tag into env '$envname'.
To use it: conda activate $envname
Remember to restart your dissect process in the new env.
EOF

if [ "$install_systemd" -ne 1 ]; then
  cat <<EOF
Tip: rerun with -s to also (re)install the wd40-checksum systemd
--user timer pointed at this env.
EOF
fi
