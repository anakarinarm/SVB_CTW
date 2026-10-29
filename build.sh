#!/usr/bin/env bash
# Written by Claude, approved by KRM.
#
# Build the MITgcm executable of the final configuration
# (code_telescoped_west64_south72_noexch2, 16 MPI ranks).
#
# Usage:
#   MITGCM_ROOT=/path/to/MITgcm ./build.sh [build_dir]
#
# MITGCM_ROOT  MITgcm checkout at commit decd05a68 (checkpoint69k-11), see README
#              (default: $HOME/Tools/MITgcm)
# OPTFILE      genmake2 build options (default: build_options/darwin_arm64_gfortran_krm)
# build_dir    default: ctw_reruns_telescoped_west64_south72/build_binary, where the
#              run folders' mitgcmuv links point.  Pass another folder to avoid
#              replacing an existing executable.
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
MITGCM_ROOT=${MITGCM_ROOT:-$HOME/Tools/MITgcm}
OPTFILE=${OPTFILE:-$HERE/build_options/darwin_arm64_gfortran_krm}
BUILD=${1:-$HERE/ctw_reruns_telescoped_west64_south72/build_binary}
PINNED=decd05a68a16c9a10b05da054f140886b217a441

if [ ! -x "$MITGCM_ROOT/tools/genmake2" ]; then
  echo "genmake2 not found under MITGCM_ROOT=$MITGCM_ROOT" >&2
  exit 1
fi
if head=$(git -C "$MITGCM_ROOT" rev-parse HEAD 2>/dev/null); then
  if [ "$head" != "$PINNED" ] &&
     ! git -C "$MITGCM_ROOT" diff --quiet "$PINNED" HEAD -- model pkg eesupp tools/genmake2 2>/dev/null; then
    echo "WARNING: MITGCM_ROOT is not at $PINNED; results may differ." >&2
  fi
fi

mkdir -p "$BUILD"
cd "$BUILD"
"$MITGCM_ROOT/tools/genmake2" -rootdir="$MITGCM_ROOT" \
  -mods="$HERE/code_telescoped_west64_south72_noexch2" -of="$OPTFILE" -mpi
make depend
make -j 8
echo "Built $BUILD/mitgcmuv"
