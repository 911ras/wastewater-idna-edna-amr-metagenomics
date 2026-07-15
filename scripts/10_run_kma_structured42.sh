#!/usr/bin/env bash

set -u
set -o pipefail

BASE="$HOME/ww_amr_batch"
DESIGN="$BASE/metadata/analysis_design.tsv"
TRIMMED="$BASE/trimmed"
DB="$BASE/kma_db/resfinder_3206"

OUTDIR="$BASE/kma_structured42"
LOGDIR="$BASE/logs/kma_structured42"

mkdir -p "$OUTDIR" "$LOGDIR"

SAMPLE_LIST="$OUTDIR/structured42_samples.txt"
STATUS="$OUTDIR/kma_batch_status.tsv"

python3 <<'PY'
from pathlib import Path
import csv

base = Path.home() / "ww_amr_batch"
design = base / "metadata" / "analysis_design.tsv"
output = base / "kma_structured42" / "structured42_samples.txt"

samples = []

with design.open(newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["INFERENTIAL_USE"] == "YES":
            samples.append(row["RUN"])

if len(samples) != 42:
    raise SystemExit(
        f"Expected 42 structured samples; found {len(samples)}"
    )

output.write_text(
    "\n".join(samples) + "\n"
)

print(f"Structured samples written: {len(samples)}")
PY

echo -e \
"SAMPLE\tSTATUS\tQUERY_FRAGMENTS\tRESULT_ROWS" \
> "$STATUS"

TOTAL=$(wc -l < "$SAMPLE_LIST")
COUNT=0

while IFS= read -r SRR; do

    COUNT=$((COUNT + 1))

    R1="$TRIMMED/${SRR}_1.trimmed.fastq.gz"
    R2="$TRIMMED/${SRR}_2.trimmed.fastq.gz"

    PREFIX="$OUTDIR/$SRR"
    LOG="$LOGDIR/${SRR}.log"

    echo
    echo "============================================================"
    echo "[$COUNT/$TOTAL] $SRR"
    echo "============================================================"

    if [ ! -s "$R1" ] || [ ! -s "$R2" ]; then

        echo "ERROR: paired FASTQ missing"

        echo -e \
        "${SRR}\tMISSING_FASTQ\tNA\tNA" \
        >> "$STATUS"

        continue
    fi

    # Resume completed results
    if [ -s "${PREFIX}.res" ] \
        && [ -s "$LOG" ] \
        && grep -q "KMA mapping done" "$LOG"
    then

        echo "EXISTS: completed KMA result"

    else

        rm -f \
        "${PREFIX}.res" \
        "${PREFIX}.aln" \
        "${PREFIX}.fsa" \
        "${PREFIX}.frag.gz"

        echo "Started: $(date)"

        /usr/bin/time -v \
        kma \
        -ipe "$R1" "$R2" \
        -t_db "$DB" \
        -o "$PREFIX" \
        -t 8 \
        > "$LOG" 2>&1

        EXIT_CODE=$?

        echo "Finished: $(date)"

        if [ "$EXIT_CODE" -ne 0 ] \
            || [ ! -s "${PREFIX}.res" ]
        then

            echo "FAILED: $SRR"

            echo -e \
            "${SRR}\tFAILED\tNA\tNA" \
            >> "$STATUS"

            continue
        fi

    fi

    QUERY_FRAGMENTS=$(
        grep \
        "Total number of query fragment after trimming" \
        "$LOG" \
        | tail -1 \
        | awk '{print $NF}'
    )

    [ -z "$QUERY_FRAGMENTS" ] \
        && QUERY_FRAGMENTS="NA"

    RESULT_ROWS=$(
        tail -n +2 "${PREFIX}.res" \
        | wc -l
    )

    echo "Query fragments: $QUERY_FRAGMENTS"
    echo "Raw KMA rows: $RESULT_ROWS"
    echo "OK: $SRR"

    echo -e \
    "${SRR}\tOK\t${QUERY_FRAGMENTS}\t${RESULT_ROWS}" \
    >> "$STATUS"

done < "$SAMPLE_LIST"

echo
echo "============================================================"
echo "KMA STRUCTURED-42 BATCH FINISHED"
echo "$(date)"
echo "============================================================"
