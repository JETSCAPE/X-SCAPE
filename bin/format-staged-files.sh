#!/usr/bin/env bash

set -euo pipefail

if ! command -v clang-format >/dev/null 2>&1 || ! command -v cmake-format >/dev/null 2>&1; then
  cat >&2 <<'EOF'
Formatter tools not found. Install them:
  Ubuntu/Debian: sudo apt install clang-format-14 cmake-format
EOF
  exit 1
fi

repo_root="$(git rev-parse --show-toplevel)"
cd "${repo_root}"

mapfile -t staged_files < <(git diff --cached --name-only --diff-filter=ACMR 2>/dev/null || true)
if [[ ${#staged_files[@]} -eq 0 ]]; then
  exit 0
fi

cpp_ext='(h|hpp|c|cc|cpp)$'
cmake_re='(^|/)CMakeLists\.txt$|\.cmake$'

selected_files=()
for f in "${staged_files[@]}"; do
  # Skip build folder and external packages
  if [[ "${f}" == build/* || "${f}" == external_packages/* ]]; then
    continue
  fi

  # Include files matching CMakeLists.txt or .cmake
  if [[ "${f}" =~ ${cmake_re} ]]; then
    selected_files+=("${f}")
    continue
  fi

  # Include source files from src/ and examples/
  if [[ "${f}" =~ \.${cpp_ext} ]] && [[ "${f}" == src/* || "${f}" == examples/* ]]; then
    selected_files+=("${f}")
  fi
done

if [[ ${#selected_files[@]} -eq 0 ]]; then
  exit 0
fi

# format the selected files
for f in "${selected_files[@]}"; do
  if [[ "${f}" =~ ${cmake_re} ]]; then
    cmake-format -i "${f}"
  else
    clang-format -i -style=file "${f}"
  fi
done

# restage if changed
if ! git diff --quiet -- "${selected_files[@]}"; then
  git add -- "${selected_files[@]}"
  cat >&2 <<'EOF'
Formatting hook changed and restaged files.
EOF
fi

exit 0
