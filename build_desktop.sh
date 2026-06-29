#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# Build a standalone desktop app for the Genome VCF Evaluator.
#
# Usage:
#   Mac:     ./build_desktop.sh          → produces  dist/Genome Analyzer.app
#   Windows: bash build_desktop.sh       → produces  dist/Genome Analyzer/Genome Analyzer.exe
#            (run in Git Bash or WSL on a Windows machine)
#
# Prerequisites:
#   pip3 install pyinstaller
#   pip3 install -r requirements.txt
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "$0")"

APP_NAME="Genome Analyzer"

# ── Ensure PyInstaller is installed ──────────────────────────────────────────
if ! python3 -m PyInstaller --version &>/dev/null; then
    echo "Installing PyInstaller…"
    pip3 install pyinstaller
fi

# ── Collect all project Python files to bundle as data ───────────────────────
DATA_FILES=(
    app.py
    chat_engine.py
    claude_analyzer.py
    database_lookup.py
    gene_lookup.py
    gene_panels.py
    main.py
    profile_manager.py
    report.py
    vcf_parser.py
)

ADD_DATA_ARGS=""
SEP=":"
# Windows uses ";" as the PyInstaller data separator
if [[ "$OSTYPE" == "msys" || "$OSTYPE" == "cygwin" || "$OSTYPE" == "win32" ]]; then
    SEP=";"
fi

for f in "${DATA_FILES[@]}"; do
    if [[ -f "$f" ]]; then
        ADD_DATA_ARGS="$ADD_DATA_ARGS --add-data ${f}${SEP}."
    fi
done

# Bundle Streamlit config directory
if [[ -d ".streamlit" ]]; then
    ADD_DATA_ARGS="$ADD_DATA_ARGS --add-data .streamlit${SEP}.streamlit"
fi

# ── Platform-specific flags ──────────────────────────────────────────────────
PLATFORM_FLAGS=""
if [[ "$OSTYPE" == "darwin"* ]]; then
    PLATFORM_FLAGS="--windowed"   # creates a .app bundle on macOS
    if [[ -f icon.icns ]]; then
        PLATFORM_FLAGS="$PLATFORM_FLAGS --icon icon.icns"
    fi
elif [[ -f icon.ico ]]; then
    PLATFORM_FLAGS="--icon icon.ico"
fi

# ── Run PyInstaller ──────────────────────────────────────────────────────────
echo ""
echo "Building ${APP_NAME}…"
echo "This may take a few minutes on first build."
echo ""

python3 -m PyInstaller \
    --name "${APP_NAME}" \
    --noconfirm \
    --clean \
    ${PLATFORM_FLAGS} \
    --collect-all streamlit \
    --collect-all cyvcf2 \
    --hidden-import anthropic \
    --hidden-import anthropic.resources \
    --hidden-import anthropic._streaming \
    --hidden-import requests \
    --hidden-import pandas \
    --hidden-import click \
    --hidden-import rich \
    --hidden-import altair \
    --hidden-import pyarrow \
    ${ADD_DATA_ARGS} \
    launcher.py

echo ""
echo "═══════════════════════════════════════════════════════════════"
if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "  ✅  Build complete: dist/${APP_NAME}.app"
    echo ""
    echo "  To run:  open \"dist/${APP_NAME}.app\""
    echo "  To distribute: zip the .app and share it."
else
    echo "  ✅  Build complete: dist/${APP_NAME}/${APP_NAME}.exe"
    echo ""
    echo "  To run:  dist\\${APP_NAME}\\${APP_NAME}.exe"
    echo "  To distribute: zip the folder and share it."
fi
echo "═══════════════════════════════════════════════════════════════"
