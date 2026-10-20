#!/usr/bin/env bash
# Build the modified PYTHIA 6.4.25 of JEWEL as a static, position-independent
# library lib/libpythia6jewel.a. Needed only where the JETSCAPE build has no
# Fortran compiler (e.g. the JETSCAPE docker image): run this with any gfortran
# (libgfortran.so.5) on the host; the result links against the libgfortran of
# the image.
#
#   ./build_pythia6.sh [gfortran]
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fc="${1:-${FC:-gfortran}}"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$here/lib"
"$fc" -O2 -fPIC -w -ffp-contract=off -c "$here/fortran/pythia6425mod-lhapdf6.f" -o "$tmp/pythia6.o"
"$fc" -O2 -fPIC -ffp-contract=off -J "$tmp" -c "$here/fortran/fio.f90" -o "$tmp/fio.o"
rm -f "$here/lib/libpythia6jewel.a"
ar rcs "$here/lib/libpythia6jewel.a" "$tmp/pythia6.o" "$tmp/fio.o"
echo "built $here/lib/libpythia6jewel.a with $("$fc" --version | head -1)"
