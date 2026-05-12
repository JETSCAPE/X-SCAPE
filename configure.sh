#!/usr/bin/env bash
# Interactive CMake configurator for X-SCAPE + js-contrib.
# Requires: dialog  (macOS: brew install dialog  |  Ubuntu: sudo apt install dialog)

set -euo pipefail

# ── Color scheme (dark, no blue) ───────────────────────────────────────────
_DIALOGRC=$(mktemp)
trap 'rm -f "$_DIALOGRC"' EXIT
cat > "$_DIALOGRC" <<'EOF'
use_shadow = ON
use_colors = ON
screen_color        = (WHITE,BLACK,OFF)
dialog_color        = (WHITE,BLACK,OFF)
title_color         = (YELLOW,BLACK,ON)
border_color        = (WHITE,BLACK,ON)
border2_color       = (WHITE,BLACK,ON)
button_active_color       = (BLACK,WHITE,ON)
button_inactive_color     = (WHITE,BLACK,OFF)
button_key_active_color   = (BLACK,YELLOW,ON)
button_key_inactive_color = (RED,BLACK,ON)
button_label_active_color = (BLACK,WHITE,ON)
button_label_inactive_color = (WHITE,BLACK,OFF)
inputbox_color        = (WHITE,BLACK,OFF)
inputbox_border_color = (YELLOW,BLACK,ON)
inputbox_border2_color = (YELLOW,BLACK,ON)
searchbox_color       = (WHITE,BLACK,OFF)
searchbox_title_color = (YELLOW,BLACK,ON)
searchbox_border_color = (YELLOW,BLACK,ON)
position_indicator_color = (YELLOW,BLACK,ON)
menubox_color         = (WHITE,BLACK,OFF)
menubox_border_color  = (YELLOW,BLACK,ON)
menubox_border2_color = (YELLOW,BLACK,ON)
item_color            = (WHITE,BLACK,OFF)
item_selected_color   = (BLACK,WHITE,ON)
tag_color             = (GREEN,BLACK,ON)
tag_selected_color    = (BLACK,GREEN,ON)
tag_key_color         = (GREEN,BLACK,ON)
tag_key_selected_color = (BLACK,GREEN,ON)
check_color           = (GREEN,BLACK,ON)
check_selected_color  = (BLACK,GREEN,ON)
uarrow_color          = (YELLOW,BLACK,ON)
darrow_color          = (YELLOW,BLACK,ON)
form_active_text_color = (BLACK,WHITE,ON)
form_text_color       = (WHITE,BLACK,OFF)
form_item_readonly_color = (WHITE,BLACK,ON)
gauge_color           = (BLACK,WHITE,ON)
EOF
export DIALOGRC="$_DIALOGRC"

# ── Locate TUI backend ──────────────────────────────────────────────────────
if command -v dialog &>/dev/null; then
    TUI=dialog
elif command -v whiptail &>/dev/null; then
    TUI=whiptail
else
    echo "ERROR: neither 'dialog' nor 'whiptail' is installed."
    echo "  macOS:  brew install dialog"
    echo "  Ubuntu: sudo apt install dialog"
    exit 1
fi

# Helper: run dialog and capture output from stdout (handles dialog's stderr default)
dlg() {
    local result
    result=$("$TUI" --stdout "$@") || { clear; echo "Cancelled."; exit 0; }
    printf '%s' "$result"
}

# ── Locate source root (directory containing this script) ──────────────────
SOURCE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Build directory ────────────────────────────────────────────────────────
BUILD_DIR=$(dlg --title "X-SCAPE Configurator" \
    --inputbox "Build directory (relative paths resolved from source root):" \
    8 65 "build")

[[ "$BUILD_DIR" = /* ]] || BUILD_DIR="${SOURCE_DIR}/${BUILD_DIR}"

# ── Build type ─────────────────────────────────────────────────────────────
BUILD_TYPE=$(dlg --title "X-SCAPE Configurator" \
    --menu "Select build type:" 12 60 4 \
    Release        "Optimised, no debug info  (recommended)" \
    Debug          "No optimisation, full debug info" \
    RelWithDebInfo "Optimised + debug info" \
    MinSizeRel     "Optimise for binary size")

# ── Main options ───────────────────────────────────────────────────────────
MAIN_CHOICES=$(dlg --title "X-SCAPE Configurator" \
    --checklist "Toggle options with SPACE, confirm with ENTER:" \
    24 76 13 \
    unittests      "Build all unit tests"                                   on  \
    USE_FREESTREAM "freestream-milne pre-equilibrium"                       off \
    USE_3DGlauber  "3D MC-Glauber initial state"                            off \
    USE_CLVISC     "CLVisc viscous hydro  (requires OpenCL)"                off \
    USE_IPGLASMA   "IP-Glasma initial state"                                off \
    USE_MUSIC      "MUSIC bulk hydrodynamics"                               off \
    USE_ISS        "iSS soft particlization  (requires MUSIC)"              off \
    USE_SMASH      "SMASH afterburner  (requires MUSIC + ISS)"              off \
    USE_ROOT       "ROOT libraries and output"                              off \
    USE_JS_CONTRIB "js-contrib extensions  (unlocks sub-options below)"    off)

# ── js-contrib sub-options (shown only when USE_JS_CONTRIB is selected) ────
JS_CHOICES=""
if echo " $MAIN_CHOICES " | grep -qw "USE_JS_CONTRIB"; then
    JS_CHOICES=$(dlg --title "js-contrib sub-options" \
        --checklist "Select js-contrib modules to build:" \
        12 76 2 \
        USE_JS_FNO_HYDRO   "FnoHydro neural-net hydro  (needs libtorch ~2 GB + ROOT)" off \
        USE_JS_PYJETSCAPE  "PyJetscape pybind11 Python bindings"                      off)
fi

# ── Check source availability and offer to run get_*.sh ───────────────────
EP="${SOURCE_DIR}/external_packages"

# option  source-dir-to-check          get script (relative to external_packages/)
declare -a PKG_OPTS PKG_DIRS PKG_SCRIPTS
PKG_OPTS=(   USE_FREESTREAM  USE_3DGlauber  USE_CLVISC  USE_IPGLASMA  USE_MUSIC  USE_ISS  USE_SMASH  USE_JS_CONTRIB )
PKG_DIRS=(   freestream-milne 3dMCGlauber  PyVisc      ipglasma      music      iSS      smash/smash_code  js-contrib )
PKG_SCRIPTS=(get_freestream-milne.sh get_3dglauber.sh get_clvisc.sh get_ipglasma.sh get_music.sh get_iSS.sh get_smash.sh get_js_contrib.sh)

MISSING_OPTS=()
MISSING_LABELS=()
for i in "${!PKG_OPTS[@]}"; do
    opt="${PKG_OPTS[$i]}"
    dir="${PKG_DIRS[$i]}"
    if echo " $MAIN_CHOICES $JS_CHOICES " | grep -qw "$opt"; then
        if [[ ! -d "${EP}/${dir}" ]]; then
            MISSING_OPTS+=("$opt")
            MISSING_LABELS+=("${opt}  →  ${PKG_SCRIPTS[$i]}")
        fi
    fi
done

if [[ ${#MISSING_OPTS[@]} -gt 0 ]]; then
    MISSING_MSG="The following enabled packages have no source yet:\n\n"
    for lbl in "${MISSING_LABELS[@]}"; do
        MISSING_MSG+="  • ${lbl}\n"
    done
    MISSING_MSG+="\nHow do you want to proceed?"

    DOWNLOAD_ACTION=$(dlg --title "Missing Sources" \
        --menu "$MISSING_MSG" 22 72 2 \
        download "Download all missing packages now" \
        skip     "Skip — I will handle downloads manually")

    if [[ "$DOWNLOAD_ACTION" == "download" ]]; then
        clear
        for i in "${!PKG_OPTS[@]}"; do
            opt="${PKG_OPTS[$i]}"
            script="${PKG_SCRIPTS[$i]}"
            if echo " ${MISSING_OPTS[*]} " | grep -qw "$opt"; then
                echo "==> Fetching ${opt}  (${script}) ..."
                (cd "${EP}" && bash "${script}")
                echo ""
            fi
        done
    else
        clear
        echo "WARNING: CMake will abort unless the missing sources are present."
        echo "Run the following from external_packages/ before building:"
        for lbl in "${MISSING_LABELS[@]}"; do
            echo "  ${lbl}"
        done
        echo ""
    fi
fi

# ── Assemble -D flags ──────────────────────────────────────────────────────
ALL_OPTIONS=(
    unittests
    USE_FREESTREAM USE_3DGlauber USE_CLVISC USE_IPGLASMA
    USE_MUSIC USE_ISS USE_SMASH USE_ROOT
    USE_JS_CONTRIB USE_JS_FNO_HYDRO USE_JS_PYJETSCAPE
)

CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=${BUILD_TYPE}"
for opt in "${ALL_OPTIONS[@]}"; do
    if echo " $MAIN_CHOICES $JS_CHOICES " | grep -qw "$opt"; then
        CMAKE_FLAGS+=" -D${opt}=ON"
    else
        CMAKE_FLAGS+=" -D${opt}=OFF"
    fi
done

FULL_CMD="cmake -S \"${SOURCE_DIR}\" -B \"${BUILD_DIR}\" ${CMAKE_FLAGS}"

# ── Confirm and run ────────────────────────────────────────────────────────
if "$TUI" --stdout --title "Confirm" \
    --yesno "Run the following command?\n\n${FULL_CMD}" \
    14 80; then
    clear
    echo "==> Build directory: ${BUILD_DIR}"
    mkdir -p "${BUILD_DIR}"
    echo "==> Running cmake..."
    echo ""
    eval "$FULL_CMD"
    echo ""
    echo "==> Done. To build:"
    echo "    cmake --build \"${BUILD_DIR}\" -j\$(nproc 2>/dev/null || sysctl -n hw.logicalcpu)"
else
    clear
    echo "Command to run manually:"
    echo ""
    echo "  mkdir -p \"${BUILD_DIR}\""
    echo "  $FULL_CMD"
    echo ""
    echo "Then build with:"
    echo "  cmake --build \"${BUILD_DIR}\" -j\$(nproc 2>/dev/null || sysctl -n hw.logicalcpu)"
fi
