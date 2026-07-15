#!/usr/bin/env bash

set -u
set -o pipefail

BASE="$HOME/ww_amr_batch"

TRIMMED="$BASE/trimmed"
DB="$BASE/kma_db/resfinder_3206"

ORIGINAL="$BASE/kma_structured42"
OUTDIR="$BASE/kma_mapstat42"
LOGDIR="$BASE/logs/kma_mapstat42"

SAMPLE_LIST="$ORIGINAL/structured42_samples.txt"
STATUS="$OUTDIR/kma_mapstat_status.tsv"

mkdir -p "$OUTDIR" "$LOGDIR"

echo -e \
"SAMPLE\tSTATUS\tMAPSTAT_FRAGMENTS\tORIGINAL_FRAGMENTS\tRES_MATCH" \
> "$STATUS"

TOTAL=$(wc -l < "$SAMPLE_LIST")
COUNT=0

while IFS= read -r SRR; do

    COUNT=$((COUNT + 1))

    R1="$TRIMMED/${SRR}_1.trimmed.fastq.gz"
    R2="$TRIMMED/${SRR}_2.trimmed.fastq.gz"

    PREFIX="$OUTDIR/$SRR"
    LOG="$LOGDIR/${SRR}.log"

    ORIGINAL_RES="$ORIGINAL/${SRR}.res"

    echo
    echo "============================================================"
    echo "[$COUNT/$TOTAL] $SRR"
    echo "============================================================"

    if [ ! -s "$R1" ] || [ ! -s "$R2" ]; then

        echo "MISSING FASTQ: $SRR"

        echo -e \
        "${SRR}\tMISSING_FASTQ\tNA\tNA\tNA" \
        >> "$STATUS"

        continue
    fi

    if [ -s "${PREFIX}.mapstat" ]; then

        echo "EXISTS: mapstat"

    else

        rm -f \
        "${PREFIX}.res" \
        "${PREFIX}.aln" \
        "${PREFIX}.fsa" \
        "${PREFIX}.frag.gz" \
        "${PREFIX}.mapstat"

        echo "Started: $(date)"

        /usr/bin/time -v \
        kma \
        -ipe "$R1" "$R2" \
        -t_db "$DB" \
        -o "$PREFIX" \
        -t 8 \
        -ef \
        > "$LOG" 2>&1

        EXIT_CODE=$?

        echo "Finished: $(date)"

        if [ "$EXIT_CODE" -ne 0 ] \
            || [ ! -s "${PREFIX}.mapstat" ] \
            || [ ! -s "${PREFIX}.res" ]
        then

            echo "FAILED: $SRR"

            echo -e \
            "${SRR}\tFAILED\tNA\tNA\tNA" \
            >> "$STATUS"

            continue
        fi

    fi

    MAPSTAT_FRAGMENTS=$(
        awk '
        $1=="##" && $2=="fragmentCount" {
            print $3
        }
        ' "${PREFIX}.mapstat" \
        | tail -1
    )

    ORIGINAL_FRAGMENTS=$(
        awk -F'\t' -v s="$SRR" '
        NR>1 && $1==s {
            print $3
        }
        ' "$ORIGINAL/kma_batch_status.tsv"
    )

    if cmp -s \
        "${PREFIX}.res" \
        "$ORIGINAL_RES"
    then
        RES_MATCH="YES"
    else
        RES_MATCH="NO"
    fi

    echo "Mapstat fragments: $MAPSTAT_FRAGMENTS"
    echo "Original fragments: $ORIGINAL_FRAGMENTS"
    echo "RES match: $RES_MATCH"

    echo -e \
    "${SRR}\tOK\t${MAPSTAT_FRAGMENTS}\t${ORIGINAL_FRAGMENTS}\t${RES_MATCH}" \
    >> "$STATUS"

    # Existing .res files are already safely stored.
    # Retain only mapstat from this new run.
    rm -f \
    "${PREFIX}.res" \
    "${PREFIX}.aln" \
    "${PREFIX}.fsa" \
    "${PREFIX}.frag.gz"

    echo "OK: $SRR"

done < "$SAMPLE_LIST"

echo
echo "============================================================"
echo "KMA MAPSTAT-42 BATCH FINISHED"
echo "$(date)"
echo "============================================================"

