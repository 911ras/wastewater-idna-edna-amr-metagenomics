#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

THREADS=10
MIN_CONTIG=500

TRIM="$HOME/ww_amr_batch/trimmed"
ASM="$HOME/ww_amr_batch/assembly"
CONTIGS="$HOME/ww_amr_batch/contigs"
AMR="$HOME/ww_amr_batch/amr_results"
LOG="$HOME/ww_amr_batch/logs"

mkdir -p "$ASM" "$CONTIGS" "$AMR" "$LOG"

cd "$TRIM"

for r1 in *_1.trimmed.fastq.gz; do
    sample="${r1%_1.trimmed.fastq.gz}"
    r2="${sample}_2.trimmed.fastq.gz"

    if [[ ! -f "$r2" ]]; then
        echo "Skipping $sample: R2 missing"
        continue
    fi

    if [[ -f "$AMR/${sample}_resfinder.tsv" ]]; then
        echo "Already completed ResFinder for $sample"
        continue
    fi

    echo "======================================"
    echo "Assembling paired-end sample: $sample"
    echo "Started: $(date)"
    echo "======================================"

    rm -rf "$ASM/${sample}_megahit"

    megahit \
      -1 "$TRIM/${sample}_1.trimmed.fastq.gz" \
      -2 "$TRIM/${sample}_2.trimmed.fastq.gz" \
      -o "$ASM/${sample}_megahit" \
      --min-contig-len "$MIN_CONTIG" \
      -t "$THREADS" \
      > "$LOG/${sample}_megahit.log" 2>&1

    contigs="$ASM/${sample}_megahit/final.contigs.fa"

    if [[ -f "$contigs" ]]; then
        cp "$contigs" "$CONTIGS/${sample}.contigs.fa"

        seqkit stats "$CONTIGS/${sample}.contigs.fa" \
          > "$AMR/${sample}_assembly_stats.txt"

        abricate --db resfinder --minid 80 --mincov 60 "$CONTIGS/${sample}.contigs.fa" \
          > "$AMR/${sample}_resfinder.tsv"

        rm -rf "$ASM/${sample}_megahit"

        echo "Completed $sample at $(date)"
    else
        echo "No contig file found for $sample" | tee -a "$LOG/failed_samples.log"
        rm -rf "$ASM/${sample}_megahit"
    fi
done

for f in *.trimmed.fastq.gz; do
    base="${f%.trimmed.fastq.gz}"

    if [[ "$base" == *_1 || "$base" == *_2 ]]; then
        continue
    fi

    sample="$base"

    if [[ -f "$AMR/${sample}_resfinder.tsv" ]]; then
        echo "Already completed ResFinder for $sample"
        continue
    fi

    echo "======================================"
    echo "Assembling single-end sample: $sample"
    echo "Started: $(date)"
    echo "======================================"

    rm -rf "$ASM/${sample}_megahit"

    megahit \
      -r "$TRIM/${sample}.trimmed.fastq.gz" \
      -o "$ASM/${sample}_megahit" \
      --min-contig-len "$MIN_CONTIG" \
      -t "$THREADS" \
      > "$LOG/${sample}_megahit.log" 2>&1

    contigs="$ASM/${sample}_megahit/final.contigs.fa"

    if [[ -f "$contigs" ]]; then
        cp "$contigs" "$CONTIGS/${sample}.contigs.fa"

        seqkit stats "$CONTIGS/${sample}.contigs.fa" \
          > "$AMR/${sample}_assembly_stats.txt"

        abricate --db resfinder --minid 80 --mincov 60 "$CONTIGS/${sample}.contigs.fa" \
          > "$AMR/${sample}_resfinder.tsv"

        rm -rf "$ASM/${sample}_megahit"

        echo "Completed $sample at $(date)"
    else
        echo "No contig file found for $sample" | tee -a "$LOG/failed_samples.log"
        rm -rf "$ASM/${sample}_megahit"
    fi
done

abricate --summary "$AMR"/*_resfinder.tsv > "$AMR/all_resfinder_summary.tsv"

echo "ResFinder-only AMR screening completed at $(date)."
