#!/usr/bin/env bash
###############################################################################
# Build the X-SCAPE documentation site locally.
#
#   MkDocs (this site + C++ API) -> ./site
#
# The C++ API reference is produced by the mkdoxy plugin, which runs Doxygen
# *during* the MkDocs build and renders native Material pages under /api/.
# There is no separate Doxygen step or standalone HTML site any more.
#
# Run from the repository root:
#     ./docs/build_docs.sh          # build into ./site
#     ./docs/build_docs.sh serve    # build, then live-preview at :8000
#
# Requirements:
#   - Python 3 + pip  (installs mkdocs-material + mkdoxy from docs/requirements.txt)
#   - doxygen + graphviz (`dot`) on PATH — required by mkdoxy for the API
#
# The GitHub Actions deploy workflow is intentionally disabled for now
# (.github/workflows/docs.yml.disabled); this script is the supported way to
# build the docs until the JETSCAPE integration is decided.
###############################################################################
set -euo pipefail

# Resolve repo root (parent of this script's directory).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${ROOT_DIR}"

echo ">> Installing MkDocs requirements ..."
python3 -m pip install -q -r docs/requirements.txt

if ! command -v doxygen >/dev/null 2>&1; then
  echo "!! doxygen not found on PATH — the mkdoxy plugin needs it for the C++ API."
  echo "   Install it, e.g. 'brew install doxygen graphviz' or"
  echo "   'sudo apt-get install doxygen graphviz', then re-run."
  exit 1
fi

echo ">> Building MkDocs site (incl. Doxygen API via mkdoxy) into ./site ..."
python3 -m mkdocs build --strict --site-dir site

echo ">> Done. Open ./site/index.html in a browser."

if [[ "${1:-}" == "serve" ]]; then
  echo ">> Starting live preview at http://localhost:8000 (Ctrl-C to stop) ..."
  python3 -m mkdocs serve
fi
