#!/bin/bash

export BLASTDB_LMDB_MAP_SIZE=0
export BLASTDB_LMDB_MAP=0

# Output directories
SPLIT_DIR="split_batches"
RESULTS_DIR="blast_results"
BEST_HITS_DIR="best_hits"
LOG_FILE="blastx_remote_log.txt"

# Create directories explicitly
mkdir -p "${RESULTS_DIR}" "${BEST_HITS_DIR}"

# Start logging
{
    echo "=============================================="
    echo " Remote BLASTX annotation started: $(date)"
    echo "=============================================="
} | tee "${LOG_FILE}"

START_TIME=$(date +%s)

# Run BLASTX remotely for each batch
for batch in "${SPLIT_DIR}"/*.fasta; do
    batch_name=$(basename "${batch}" .fasta)
    output_file="${RESULTS_DIR}/${batch_name}_blastx_results.txt"

    # Skip if already done
    if [[ -s "${output_file}" ]]; then
        echo "[⏭️ SKIP] Already completed: ${batch_name}" | tee -a "${LOG_FILE}"
        continue
    fi

    # Confirm file paths explicitly (for debugging)
    echo "[DEBUG] Batch: ${batch}, Output: ${output_file}" | tee -a "${LOG_FILE}"

    blastx -query "${batch}" \
           -db /media/eugenia/Main/Eugenia_projects/ncbi_taxonomy/viridiplantae_db/nr_viridiplantae_db \
           -evalue 1e-10 \
           -num_threads 28 \
           -max_target_seqs 5 \
           -outfmt "6 qseqid sseqid evalue bitscore stitle" \
           -out "${output_file}"

    # Check BLASTX status and output file size explicitly
    if [[ $? -eq 0 && -s "${output_file}" ]]; then
        echo "[✅ SUCCESS] BLASTX completed successfully for batch: ${batch_name}" | tee -a "${LOG_FILE}"
    else
        echo "[❌ ERROR] BLASTX failed or output empty for batch: ${batch_name}" | tee -a "${LOG_FILE}"
        rm -f "${output_file}"  # Remove empty or failed results
        continue
    fi

    # Filter for best hits per transcript (sorted by bitscore)
    best_hit_output="${BEST_HITS_DIR}/${batch_name}_best_hits.txt"
    sort -k1,1 -k4,4nr "${output_file}" | awk '!seen[$1]++' > "${best_hit_output}"

    echo "[INFO] Best hits filtered for batch: ${batch_name}" | tee -a "${LOG_FILE}"

done

# Combine results afterwards:
cat "${BEST_HITS_DIR}"/*_best_hits.txt > "${BEST_HITS_DIR}/nr_best_hits.txt"

# End timing
END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))

# Log total runtime
{
    echo "=============================================="
    echo " Remote BLASTX annotation finished: $(date)"
    echo " Total elapsed time: ${ELAPSED_TIME} seconds (~$((ELAPSED_TIME / 60)) minutes)"
    echo " Final results saved at: ${BEST_HITS_DIR}/nr_best_hits.txt"
    echo "=============================================="
} | tee -a "${LOG_FILE}"
