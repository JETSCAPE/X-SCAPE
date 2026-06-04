#!/usr/bin/env bash
###############################################################################
# Build the X-SCAPE documentation site locally.
#
#   MkDocs (this site)  ->  ./site
#   Doxygen (C++ API)   ->  ./site/api/cpp
#
# Run from the repository root:
#     ./docs/build_docs.sh          # build into ./site
#     ./docs/build_docs.sh serve    # build, then live-preview at :8000
#
# Requirements:
#   - Python 3 + pip  (installs mkdocs-material from docs/requirements.txt)
#   - doxygen + graphviz (`dot`) for the API reference (optional; skipped if
#     doxygen is not on PATH)
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

echo ">> Building MkDocs site into ./site ..."
python3 -m mkdocs build --strict --site-dir site

if command -v doxygen >/dev/null 2>&1; then
  echo ">> Building Doxygen API reference into ./site/api/cpp ..."
  doxygen docs/Doxyfile
else
  echo "!! doxygen not found on PATH — skipping the C++ API reference."
  echo "   Install it (e.g. 'sudo apt-get install doxygen graphviz') to include it."
fi

echo ">> Done. Open ./site/index.html in a browser."

if [[ "${1:-}" == "serve" ]]; then
  echo ">> Starting live preview at http://localhost:8000 (Ctrl-C to stop) ..."
  python3 -m mkdocs serve
fi
