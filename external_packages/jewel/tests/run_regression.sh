#!/bin/bash
# Regression test: the C++ port must reproduce the HepMC output of the original
# Fortran JEWEL 2.6.0 (gfortran -O2 -ffp-contract=off, macOS arm64) byte for byte.
# usage: run_regression.sh <jewel binary> <sample: pp|aar> <tests dir>
#
# expected.sha256 holds the SHA-256 of test.hepmc written by the ORIGINAL Fortran
# (built with -O2 -ffp-contract=off) for the same parameter files, with the lookup
# tables computed in the run (no table files present). To regenerate on another
# platform, for s in pp (medium-vac) and aar (medium-simple):
#   mkdir t && cp tests/$s.params.dat tests/medium.params.dat t/ && cd t
#   <reference>/jewel-2.6.0-<vac|simple> $s.params.dat && shasum -a 256 test.hepmc
set -eu
BIN=$1; S=$2; T=$3
D=$(mktemp -d); trap 'rm -rf "$D"' EXIT
cp "$T/$S.params.dat" "$T/medium.params.dat" "$D/"
cd "$D"
"$BIN" "$S.params.dat" > stdout.txt 2>&1
got=$(shasum -a 256 test.hepmc | cut -d' ' -f1)
want=$(awk -v s="$S" '$1==s{print $2}' "$T/expected.sha256")
echo "$S: got $got"
echo "$S: want $want"
[ "$got" = "$want" ]
