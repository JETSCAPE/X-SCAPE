#!/usr/bin/env bash

set -euo pipefail

if [[ ! -d .git ]]; then
  echo "Error: Run this script from the repository's root directory." >&2
  exit 1
fi

if [[ ! -x .githooks/pre-commit ]]; then
  echo "Error: Expected hook .githooks/pre-commit not found." >&2
  exit 1
fi

git config core.hooksPath .githooks
echo "Installed optional formatting hook."

if ! command -v clang-format >/dev/null 2>&1 || ! command -v cmake-format >/dev/null 2>&1; then
  cat >&2 <<'EOF'
Warning: Formatter utilities not found. Install them with your package manager:
  Ubuntu / Debian:
    sudo apt install clang-format-14 cmake-format
    Package manager commands for other distributions may differ.
  After installing the formatter utilities, re-run this setup script.
EOF
  exit 1
fi

echo "clang-format and cmake-format are installed."

