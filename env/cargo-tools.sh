#!/bin/bash
# Source-built contestants. Run once on a login node with cargo available
# (`module load rust` or rustup). Binaries land in $SB_BIN.
set -euo pipefail
source "$(dirname "$0")/../scripts/hpc/config.sh"
mkdir -p "$SB_BIN" "$SB_SRC"

# deacon-syncmer: upstream deacon plus a --scheme flag on `index build`.
# The patch is env/deacon-syncmer.patch; DEACON_FORK_URL may point at your fork
# instead once the patch is committed there.
if [ ! -x "$SB_BIN/deacon-syncmer" ]; then
  if [ ! -d "$SB_SRC/deacon" ]; then
    git clone --depth 50 "${DEACON_UPSTREAM_URL}" "$SB_SRC/deacon"
  fi
  ( cd "$SB_SRC/deacon"
    git checkout -q "${DEACON_UPSTREAM_TAG}"
    git apply --check "$REPO_DIR/env/deacon-syncmer.patch"
    git apply "$REPO_DIR/env/deacon-syncmer.patch"
    cargo build --release
    cp target/release/deacon "$SB_BIN/deacon-syncmer" )
fi

# 2b-RAD tool rows: Syn2bANI (T1) and Fast2bRAD-M (T2/T3), both Rust.
build_rust() {   # <name> <url> <tag> <binary-name-in-target>
  local name="$1" url="$2" tag="$3" bin="$4"
  [ -x "$SB_BIN/$bin" ] && return 0
  [ -d "$SB_SRC/$name" ] || git clone "$url" "$SB_SRC/$name"
  ( cd "$SB_SRC/$name" && git fetch -q && git checkout -q "$tag" && cargo build --release && cp "target/release/$bin" "$SB_BIN/$bin" ) \
    || echo "WARN: $name build failed; its rows will be skipped"
}
build_rust syn2bani   "$SYN2BANI_URL"  "$SYN2BANI_TAG"  syn2bani
build_rust fast2brad  "$FAST2BRAD_URL" "$FAST2BRAD_TAG" fast2bRAD-M

# gsearch (ProbMinHash + HNSW). Optional; the ProbMinHash row also runs via
# Dashing 2's --pminhash mode, so a failed build here loses nothing.
if [ ! -x "$SB_BIN/gsearch" ]; then
  cargo install --root "$SB_BIN/.." gsearch || echo "WARN: gsearch build failed; ProbMinHash row will use dashing2 --pminhash"
fi

# dashing2: static release binary from GitHub.
if [ ! -x "$SB_BIN/dashing2" ]; then
  curl -L -o "$SB_BIN/dashing2" "${DASHING2_URL}" && chmod +x "$SB_BIN/dashing2" \
    || echo "WARN: dashing2 download failed; set DASHING2_URL to a release asset for x86_64 linux"
fi

echo "installed:"; ls -la "$SB_BIN"
