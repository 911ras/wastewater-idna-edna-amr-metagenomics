#!/usr/bin/env bash
set -euo pipefail

RAW="/mnt/c/sra_project/ww_amr_sra/fastq"
OUT="$HOME/ww_amr_batch/trimmed"
REPORT="$HOME/ww_amr_batch/fastp_report"
LOG="$HOME/ww_amr_batch/logs"

mkdir -p "$OUT" "$REPORT" "$LOG"

cd "$RAW"

for f in *.fastq.gz; do
    base="${f%.fastq.gz}"

    # Skip R2 files; they are handled with R1
    if [[ "$base" == *_2 ]]; then
        continue
    fi

    # Case 1: paired-end sample
    if [[ "$base" == *_1 && -f "${base%_1}_2.fastq.gz" ]]; then
        sample="${base%_1}"

        if [[ -f "$OUT/${sample}_1.trimmed.fastq.gz" && -f "$OUT/${sample}_2.trimmed.fastq.gz" ]]; then
            echo "Already done, skipping paired-end $sample"
            continue
        fi

        echo "Running paired-end fastp for $sample"

        fastp \
          -i "$RAW/${sample}_1.fastq.gz" \
          -I "$RAW/${sample}_2.fastq.gz" \
          -o "$OUT/${sample}_1.trimmed.fastq.gz" \
          -O "$OUT/${sample}_2.trimmed.fastq.gz" \
          -h "$REPORT/${sample}_fastp.html" \
          -j "$REPORT/${sample}_fastp.json" \
          --detect_adapter_for_pe \
          --thread 4 \
          > "$LOG/${sample}_fastp.log" 2>&1

    # Case 2: single-end file named SRRxxxx_1.fastq.gz but no _2 exists
    elif [[ "$base" == *_1 ]]; then
        sample="${base%_1}"

        if [[ -f "$OUT/${sample}.trimmed.fastq.gz" ]]; then
            echo "Already done, skipping single-end $sample"
            continue
        fi

        echo "Running single-end fastp for $sample from ${base}.fastq.gz"

        fastp \
          -i "$RAW/${base}.fastq.gz" \
          -o "$OUT/${sample}.trimmed.fastq.gz" \
          -h "$REPORT/${sample}_fastp.html" \
          -j "$REPORT/${sample}_fastp.json" \
          --thread 4 \
          > "$LOG/${sample}_fastp.log" 2>&1

    # Case 3: normal single-end file named SRRxxxx.fastq.gz
    else
        sample="$base"

        if [[ -f "$OUT/${sample}.trimmed.fastq.gz" ]]; then
            echo "Already done, skipping single-end $sample"
            continue
        fi

        echo "Running single-end fastp for $sample"

        fastp \
          -i "$RAW/${sample}.fastq.gz" \
          -o "$OUT/${sample}.trimmed.fastq.gz" \
          -h "$REPORT/${sample}_fastp.html" \
          -j "$REPORT/${sample}_fastp.json" \
          --thread 4 \
          > "$LOG/${sample}_fastp.log" 2>&1
    fi
done

echo "fastp trimming completed."
