#!/usr/bin/env bash
# vendor_piglet.sh — vendor the PIgLET R package for Docker builds
#
# PIgLET (https://bitbucket.org/yaarilab/piglet) is installed into the
# Docker image from vendor/piglet-built/ — a copy of the compiled package
# from a local R library. This avoids Bitbucket auth in CI and compilation
# issues in the Docker build environment.
#
# Usage:
#   bash vendor_piglet.sh              # auto-detect local piglet
#   bash vendor_piglet.sh /path/to/R/library/piglet
#
# After running:
#   git add vendor/piglet-built
#   git commit -m "vendor pre-built piglet <version>"
#   git push

set -euo pipefail

DEST="vendor/piglet-built"

# Find piglet installation
if [[ -n "${1:-}" && -d "$1" ]]; then
  SRC="$1"
else
  # Auto-detect from R
  SRC=$(Rscript -e "cat(system.file(package='piglet'))" 2>/dev/null)
fi

[[ -z "$SRC" || ! -d "$SRC" ]] && {
  echo "ERROR: piglet not found. Install it first:"
  echo "  library(devtools); install_bitbucket('yaarilab/piglet')"
  exit 1
}

VERSION=$(Rscript -e "cat(as.character(packageVersion('piglet')))" 2>/dev/null)
RVERSION=$(Rscript -e "cat(R.version\$major, '.', R.version\$minor, sep='')" 2>/dev/null)
ARCH=$(uname -m)
OS=$(uname -s)

echo "Vendoring PIgLET ${VERSION} from: ${SRC}"
echo "  Built with: R ${RVERSION} on ${OS} ${ARCH}"
echo ""

if [[ "${OS}" != "Linux" ]]; then
  echo "WARNING: You are on ${OS}. The compiled .so files will only work"
  echo "  in a Docker container if the architecture and glibc version match."
  echo "  The Docker base (rocker/r-ver:4.4.1) uses Ubuntu 22.04 x86_64."
  echo "  Build on a Linux x86_64 machine for guaranteed compatibility."
  echo ""
fi

# Check GLIBCXX requirement of the compiled .so
SO="${SRC}/libs/piglet.so"
if [[ -f "$SO" ]]; then
  GLIBCXX_MAX=$(objdump -p "$SO" 2>/dev/null | grep GLIBCXX | \
    awk '{print $NF}' | sort -V | tail -1 || echo "unknown")
  echo "  Requires GLIBCXX: ${GLIBCXX_MAX}"
  echo "  Docker base provides up to: GLIBCXX_3.4.30 (Ubuntu 22.04)"
  if [[ "$GLIBCXX_MAX" > "GLIBCXX_3.4.30" ]]; then
    echo ""
    echo "WARNING: ${GLIBCXX_MAX} > GLIBCXX_3.4.30"
    echo "  This .so may fail to load in the Docker container."
    echo "  Consider switching the Docker base to Ubuntu 24.04:"
    echo "    FROM rocker/r-ver:4.4.1 (Ubuntu 24.04 based)"
    echo "  or compile PIgLET inside the container instead."
    echo ""
  fi
fi

rm -rf "$DEST"
cp -r "$SRC" "$DEST"
echo "Vendored to: $DEST"
echo ""
echo "Next steps:"
echo "  git add vendor/piglet-built"
echo "  git commit -m 'vendor pre-built piglet ${VERSION}'"
echo "  git push"
