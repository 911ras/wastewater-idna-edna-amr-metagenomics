#!/usr/bin/env bash
set -euo pipefail

BASE="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
TMP="$BASE/MANIFEST_SHA256.tsv.tmp"

(
  cd "$BASE"

  find . \
    -path './.git' -prune -o \
    -type f \
    ! -name 'MANIFEST_SHA256.tsv' \
    ! -name 'MANIFEST_SHA256.tsv.tmp' \
    -print0 \
    | sort -z \
    | while IFS= read -r -d '' file; do
        hash="$(sha256sum "$file" | awk '{print $1}')"
        size="$(stat -c '%s' "$file")"
        printf '%s\t%s\t%s\n' \
          "$hash" "$size" "${file#./}"
      done
) > "$TMP"

mv "$TMP" "$BASE/MANIFEST_SHA256.tsv"

echo "Manifest regenerated:"
echo "$BASE/MANIFEST_SHA256.tsv"
echo "Entries: $(wc -l < "$BASE/MANIFEST_SHA256.tsv")"
