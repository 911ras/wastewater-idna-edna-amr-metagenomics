#!/usr/bin/env bash
set -euo pipefail

BASE="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
REF_DIR="$BASE/references"
mkdir -p "$REF_DIR"

DB_ROOT="${CONDA_PREFIX:-}/db/resfinder"
SOURCE="$DB_ROOT/sequences"
DEST="$REF_DIR/resfinder_3206_2026-06-29.fasta"

if [[ -z "${CONDA_PREFIX:-}" || ! -f "$SOURCE" ]]; then
    echo "ResFinder FASTA not found at: $SOURCE" >&2
    echo "Activate the exact AMR environment and rerun this script." >&2
    exit 1
fi

cp -p "$SOURCE" "$DEST"
sha256sum "$DEST" > "$REF_DIR/resfinder_3206_2026-06-29.fasta.sha256"

{
    echo "Captured from: $SOURCE"
    echo "ABRicate database record:"
    abricate --list 2>/dev/null | grep -i '^resfinder' || true
    echo "FASTA records: $(grep -c '^>' "$DEST")"
    echo "SHA256: $(sha256sum "$DEST" | awk '{print $1}')"
} > "$REF_DIR/resfinder_reference_snapshot_info.txt"

echo "Reference snapshot captured: $DEST"
echo "Records: $(grep -c '^>' "$DEST")"
echo "SHA256 file: $REF_DIR/resfinder_3206_2026-06-29.fasta.sha256"
