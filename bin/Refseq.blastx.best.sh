#!/bin/bash

# Transcriptome Annotation using BLASTX  and RefSeq DB (Best hit per transcript)

SPLIT_DIR="split_batches" #split file with into batches, too big
RESULTS_DIR="diamond_results"
BEST_HITS_DIR="best_hits"
DB_NAME="refseq_proteins.dmnd"  # DIAMOND DB
FASTA_DB="all_refseq_proteins.faa"  # Combined RefSeq protein fasta
LOG_FILE="diamond_annotation_log.txt"

# === PREPARE ===
mkdir -p "${RESULTS_DIR}" "${BEST_HITS_DIR}"

#LOG
{
    echo "=============================================="
    echo " Local DIAMOND annotation started: $(date)"
    echo "=============================================="
} | tee "${LOG_FILE}"

START_TIME=$(date +%s)

# CHECK IF DB EXISTS
if [[ ! -f "${DB_NAME}" ]]; then
    echo "[INFO] Building DIAMOND database from: ${FASTA_DB}" | tee -a "${LOG_FILE}"
    diamond makedb --in "${FASTA_DB}" -d refseq_proteins
    if [[ $? -ne 0 ]]; then
        echo "[❌ ERROR] Failed to create DIAMOND database!" | tee -a "${LOG_FILE}"
        exit 1
    fi
fi

#  PROCESS EACH BATCH
for batch in "${SPLIT_DIR}"/*.fasta; do
    batch_name=$(basename "${batch}" .fasta)
    output_file="${RESULTS_DIR}/${batch_name}_diamond.tsv"

    echo "[DEBUG] Processing: ${batch_name}" | tee -a "${LOG_FILE}"

    diamond blastx \
        -q "${batch}" \
        -d "${DB_NAME}" \
        -o "${output_file}" \
        -f 6 qseqid sseqid evalue bitscore stitle \
        --evalue 1e-10 \
        --max-target-seqs 5 \
        --threads 30

    # Check status
    if [[ $? -eq 0 ]]; then
        echo "[✅ SUCCESS] DIAMOND completed for batch: ${batch_name}" | tee -a "${LOG_FILE}"
    else
        echo "[❌ ERROR] DIAMOND failed for batch: ${batch_name}" | tee -a "${LOG_FILE}"
        continue
    fi

    # BEST HIT FILTER
    best_hit_output="${BEST_HITS_DIR}/${batch_name}_best_hits.txt"
    sort -k1,1 -k4,4nr "${output_file}" | awk '!seen[$1]++' > "${best_hit_output}"

    echo "[INFO] Best hits saved for batch: ${batch_name}" | tee -a "${LOG_FILE}"
done

# === COMBINE RESULTS ===
cat "${BEST_HITS_DIR}"/*_best_hits.txt > "${BEST_HITS_DIR}/refseq_best_hits.txt"

END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))

# === FINAL LOG ===
{
    echo "=============================================="
    echo " DIAMOND annotation finished: $(date)"
    echo " Total elapsed time: ${ELAPSED_TIME} seconds (~$((ELAPSED_TIME / 60)) minutes)"
    echo " Final results at: ${BEST_HITS_DIR}/refseq_best_hits.txt"
    echo "=============================================="
} | tee -a "${LOG_FILE}"
