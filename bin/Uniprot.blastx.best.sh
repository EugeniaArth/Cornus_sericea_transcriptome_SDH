#!/bin/bash


# Transcriptome Annotation using BLASTX against UniProt database (Best hit per transcript)

# Paths
OUTPUT_DIR="Annotation/Uniprot"
OUTPUT_FILE_RAW="${OUTPUT_DIR}/blastx_annotations_raw.txt"
OUTPUT_FILE_BEST="${OUTPUT_DIR}/blastx_annotations_best_hits.txt"
LOG_FILE="${OUTPUT_DIR}/blastx_log.txt"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Logging
{
    echo "==============================="
    echo "BLASTX Annotation Started: $(date)"
    echo "Input File: final.clust_transcripts_longest_iso"
    echo "Database: blastx/plant_db_uniprot/uniprot_plants_db"
    echo "==============================="
} | tee "$LOG_FILE"

# Run BLASTX
START_TIME=$(date +%s)

blastx -query final.clust_transcripts_longest_iso.fasta \
       -db plant_db_uniprot/uniprot_plants_db \
       -evalue 1e-10 \
       -num_threads 30 \
       -outfmt "6 qseqid sseqid evalue bitscore stitle" \
       -max_target_seqs 5 \
       -out blastx_annotations_raw.txt 2>> "$LOG_FILE"

if [[ $? -eq 0 ]]; then
    echo "BLASTX completed successfully!" | tee -a "$LOG_FILE"
else
    echo "Error: BLASTX failed. Check log: $LOG_FILE" | tee -a "$LOG_FILE"
    exit 1
fi

# Filter best hits per transcript
awk 'BEGIN {FS="\t"; OFS="\t"}
     !seen[$1]++ {print $0}' <(sort -k1,1 -k3,3g -k4,4gr blastx_annotations_raw.txt) > uniprot_best_hits.txt

if [[ $? -eq 0 ]]; then
    echo "Successfully filtered best hits." | tee -a "$LOG_FILE"
else
    echo "Error: Filtering best hits failed." | tee -a "$LOG_FILE"
    exit 1
fi

# End timing
END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))

# Logging
{
    echo "==============================="
    echo "BLASTX Finished: $(date)"
    echo "Total Time: $ELAPSED_TIME seconds (~$(($ELAPSED_TIME / 60)) minutes)"
    echo "Raw Results: blastx_annotations_raw.txt"
    echo "Best hits Results: blastx_annotations_best_hits.txt"
    echo "==============================="
} | tee -a "$LOG_FILE"
