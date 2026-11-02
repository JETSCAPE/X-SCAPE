#!/usr/bin/env bash
#
# Run a single X-SCAPE (runJetscape) instance in its own working directory,
# reading read-only assets from a shared prefix. This lets many instances run
# in parallel from ONE build/install tree, instead of one full build directory
# per concurrent run.
#
# The shared asset prefix ("data dir") is the directory that holds EOS/,
# iSS_tables/, iSS_parameters.dat, LBT-tables/, tables/, eps09/, LHAPDF_Lib/,
# nucleusConfigs/, data_table/, mcglauber.input, music_input and config/. For an
# in-build run this is simply the build directory.
#
# Read-only assets are located via:
#   * XSCAPE_DATA_DIR  -> mcglauber.input (resolved in MCGlauberWrapper)
#   * HYDROPROGRAMPATH -> MUSIC EOS tables (resolved in MUSIC)
#   * LBT_TABLES_PATH  -> LBT tables (resolved in LBT and in Matter's recoil)
#   * per-run copies of jetscape_main.xml and of the user XML with every '../'
#     path made absolute (run_xml.py, next to this script): relative to DATA_DIR,
#     else found in DATA_DIR by its trailing parts (an installed tree keeps
#     iSS_tables, iSS_parameters.dat and config/ there); the user's file is not
#     changed
#   * symlinks for the few dirs the vendored code still reads CWD-relative
# Before running, run_xml.py stops (exit 3) when the user XML has a tag the main
# XML lacks, or when an input the run uses has no file behind it; -f runs anyway.
# Writable per-run state (music_input and all simulation output) lives in the
# run directory, so concurrent runs never collide.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: run_in_workdir.sh -w RUN_DIR -u USER_XML [-d DATA_DIR] [-b RUNJETSCAPE] [-m MAIN_XML] [-f]

  -w RUN_DIR      Per-run working directory (created if missing). REQUIRED.
  -u USER_XML     User XML config. REQUIRED.
  -d DATA_DIR     Shared asset prefix. Default: $XSCAPE_DATA_DIR if set,
                  else the directory containing the runJetscape executable.
  -b RUNJETSCAPE  Path to the runJetscape executable.
                  Default: DATA_DIR/runJetscape.
  -m MAIN_XML     Main XML config. Default: DATA_DIR/config/jetscape_main.xml.
  -f              Run even if a check fails: a user XML tag the main XML lacks,
                  or an input of the run with no file behind it.

Both XMLs are copied into RUN_DIR (jetscape_main.xml, <user XML stem>.run.xml)
with their '../' paths made absolute, so any XML of the source tree or the
installed config/ runs as it is. Give inputs of your own as absolute paths.

Example (two parallel runs against one build dir):
  ./run_in_workdir.sh -d /path/to/build -w /tmp/run1 -u /path/to/build/config/jetscape_user_3DGlauber_MUSIC_iSS_SMASH_test.xml &
  ./run_in_workdir.sh -d /path/to/build -w /tmp/run2 -u /path/to/build/config/jetscape_user_3DGlauber_MUSIC_iSS_SMASH_test.xml &
  wait
EOF
}

# Portable absolute-path helpers (macOS lacks GNU realpath by default).
abspath_dir() { (cd "$1" >/dev/null 2>&1 && pwd); }
abspath_file() {
  local d b
  d=$(cd "$(dirname "$1")" >/dev/null 2>&1 && pwd) || return 1
  b=$(basename "$1")
  printf '%s/%s\n' "$d" "$b"
}

RUN_DIR=""
USER_XML=""
DATA_DIR="${XSCAPE_DATA_DIR:-}"
RUNJETSCAPE=""
MAIN_XML=""
FORCE=0

while getopts ":w:u:d:b:m:fh" opt; do
  case "$opt" in
    w) RUN_DIR="$OPTARG" ;;
    u) USER_XML="$OPTARG" ;;
    d) DATA_DIR="$OPTARG" ;;
    b) RUNJETSCAPE="$OPTARG" ;;
    m) MAIN_XML="$OPTARG" ;;
    f) FORCE=1 ;;
    h) usage; exit 0 ;;
    \?) echo "Unknown option: -$OPTARG" >&2; usage; exit 2 ;;
    :) echo "Option -$OPTARG requires an argument" >&2; usage; exit 2 ;;
  esac
done

[ -n "$RUN_DIR" ] || { echo "ERROR: -w RUN_DIR is required" >&2; usage; exit 2; }
[ -n "$USER_XML" ] || { echo "ERROR: -u USER_XML is required" >&2; usage; exit 2; }

# Resolve DATA_DIR: explicit/-d/env, else infer from the runJetscape location.
if [ -z "$DATA_DIR" ] && [ -n "$RUNJETSCAPE" ]; then
  DATA_DIR=$(dirname "$RUNJETSCAPE")
fi
[ -n "$DATA_DIR" ] || { echo "ERROR: data dir unknown; pass -d or set XSCAPE_DATA_DIR" >&2; exit 2; }
DATA_DIR=$(abspath_dir "$DATA_DIR") || { echo "ERROR: data dir not found: $DATA_DIR" >&2; exit 1; }

[ -n "$RUNJETSCAPE" ] || RUNJETSCAPE="$DATA_DIR/runJetscape"
[ -x "$RUNJETSCAPE" ] || { echo "ERROR: runJetscape not executable: $RUNJETSCAPE" >&2; exit 1; }
RUNJETSCAPE=$(abspath_file "$RUNJETSCAPE")

# Default main XML: installed layout (DATA_DIR/config) or in-build layout
# (config lives at the source root, one level up from the build dir).
if [ -z "$MAIN_XML" ]; then
  for cand in "$DATA_DIR/config/jetscape_main.xml" "$DATA_DIR/../config/jetscape_main.xml"; do
    if [ -f "$cand" ]; then MAIN_XML="$cand"; break; fi
  done
fi
[ -n "$MAIN_XML" ] && [ -f "$MAIN_XML" ] || { echo "ERROR: main XML not found (pass -m); looked in $DATA_DIR/config and $DATA_DIR/../config" >&2; exit 1; }
MAIN_XML=$(abspath_file "$MAIN_XML")

[ -f "$USER_XML" ] || { echo "ERROR: user XML not found: $USER_XML" >&2; exit 1; }
USER_XML=$(abspath_file "$USER_XML")

mkdir -p "$RUN_DIR"
RUN_DIR=$(abspath_dir "$RUN_DIR")

# Per-run copies of both XMLs with their '../' paths absolute, so the run works
# from any CWD (the shipped configs use ../external_packages/iSS/... and other
# paths relative to the build directory), and the checks (run_xml.py).
RUN_XML_PY="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_xml.py"
PYTHON=${PYTHON:-python3}
if [ -f "$RUN_XML_PY" ] && command -v "$PYTHON" >/dev/null 2>&1; then
  check=()
  [ "$FORCE" = 1 ] && check=(--no-check)
  xmls=$("$PYTHON" "$RUN_XML_PY" -d "$DATA_DIR" -m "$MAIN_XML" -u "$USER_XML" -w "$RUN_DIR" ${check[@]+"${check[@]}"}) || exit $?
  RUN_MAIN_XML=$(printf '%s\n' "$xmls" | sed -n 1p)
  RUN_USER_XML=$(printf '%s\n' "$xmls" | sed -n 2p)
else
  # without Python: only the main XML's iSS paths, as before
  echo "[run_in_workdir] WARNING: $PYTHON or $RUN_XML_PY not found: only the main XML's iSS paths are made absolute, nothing is checked" >&2
  RUN_MAIN_XML="$RUN_DIR/jetscape_main.xml"
  RUN_USER_XML="$USER_XML"
  sed -E \
    -e "s#<iSS_input_file>[^<]*</iSS_input_file>#<iSS_input_file>${DATA_DIR}/iSS_parameters.dat</iSS_input_file>#" \
    -e "s#<iSS_table_path>[^<]*</iSS_table_path>#<iSS_table_path>${DATA_DIR}/iSS_tables</iSS_table_path>#" \
    -e "s#<iSS_particle_table_path>[^<]*</iSS_particle_table_path>#<iSS_particle_table_path>${DATA_DIR}/iSS_tables</iSS_particle_table_path>#" \
    "$MAIN_XML" > "$RUN_MAIN_XML"
fi

# Per-run writable copy of the MUSIC input file (MUSIC rewrites it in place).
# The shipped config references it CWD-relative as "music_input".
if [ -f "$DATA_DIR/music_input" ]; then
  cp -f "$DATA_DIR/music_input" "$RUN_DIR/music_input"
fi

# Symlink the dirs that the vendored 3dMCGlauber / trento code still reads
# CWD-relative. These are read-only, so sharing them across runs is safe.
# LBT-tables: Matter's recoil now honours LBT_TABLES_PATH too, but binaries
# built before that read ./LBT-tables, so the link keeps them working.
for asset in tables eps09 LHAPDF_Lib nucleusConfigs data_table LBT-tables; do
  if [ -e "$DATA_DIR/$asset" ] && [ ! -e "$RUN_DIR/$asset" ]; then
    ln -s "$DATA_DIR/$asset" "$RUN_DIR/$asset"
  fi
done

export XSCAPE_DATA_DIR="$DATA_DIR"
export HYDROPROGRAMPATH="$DATA_DIR"
export LBT_TABLES_PATH="$DATA_DIR/LBT-tables"

echo "[run_in_workdir] DATA_DIR   = $DATA_DIR"
echo "[run_in_workdir] RUN_DIR    = $RUN_DIR"
echo "[run_in_workdir] runJetscape= $RUNJETSCAPE"
echo "[run_in_workdir] user XML   = $RUN_USER_XML (from $USER_XML)"
echo "[run_in_workdir] main XML   = $RUN_MAIN_XML"

cd "$RUN_DIR"
exec "$RUNJETSCAPE" "$RUN_USER_XML" "$RUN_MAIN_XML"
