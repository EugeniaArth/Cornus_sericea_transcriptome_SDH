import pandas as pd
from goatools.obo_parser import GODag

# Input / output files
input_file = "Files/eggnog_annotations.tsv"
obo_file = "db/go-basic.obo"

output_file_all = "GO_terms_with_ontology.csv"
output_file_filtered = "GO_terms_highest_level_per_transcript.csv"

# Column indices for eggNOG-mapper annotation file (no header)
TRANSCRIPT_COL = 0
GO_COL = 9  # based on sample input

print("[INFO] Reading eggNOG annotation file...")
eggnog_df = pd.read_csv(
    input_file,
    sep="\t",
    header=None,
    comment="#",
    dtype=str,
    low_memory=False
)

print(f"[INFO] Loaded {len(eggnog_df)} rows.")

print("\n[DEBUG] First 5 values in GO column:")
print(eggnog_df.iloc[:5, GO_COL])

# Extract transcript–GO pairs
go_data = []
skipped_transcripts = 0

for _, row in eggnog_df.iterrows():
    transcript = row.iloc[TRANSCRIPT_COL]

    if len(row) <= GO_COL:
        skipped_transcripts += 1
        continue

    go_field = row.iloc[GO_COL]

    if pd.isna(go_field) or str(go_field).strip() in {"", "-"}:
        skipped_transcripts += 1
        continue

    go_terms = [
        term.strip()
        for term in str(go_field).split(",")
        if term.strip().startswith("GO:")
    ]

    if not go_terms:
        skipped_transcripts += 1
        continue

    for go_id in go_terms:
        go_data.append((transcript, go_id))

print(f"\n[INFO] Parsed {len(go_data)} transcript–GO pairs.")
print(f"[INFO] Skipped {skipped_transcripts} transcripts with no GO terms.")

# Create DataFrame
go_df = pd.DataFrame(go_data, columns=["Transcript_ID", "GO_term"])

if go_df.empty:
    print("[WARNING] No GO terms were found. No output files written.")
    raise SystemExit

# Remove duplicate transcript–GO pairs if present
go_df = go_df.drop_duplicates()

print("\n[DEBUG] First few transcript–GO pairs:")
print(go_df.head())

# Load GO ontology DAG
print("\n[INFO] Loading GO DAG...")
go_dag = GODag(obo_file)

# Annotation helper functions
def get_namespace(go_id):
    term = go_dag.get(go_id)
    return term.namespace if term else "Unknown"

def get_name(go_id):
    term = go_dag.get(go_id)
    return term.name if term else "Unknown"

def get_level(go_id):
    term = go_dag.get(go_id)
    return term.level if term else -1  # use -1 for unknown terms

# Add GO metadata
print("[INFO] Annotating GO terms...")
go_df["Ontology"] = go_df["GO_term"].apply(get_namespace)
go_df["GO_name"] = go_df["GO_term"].apply(get_name)
go_df["GO_level"] = go_df["GO_term"].apply(get_level)

# Save full annotated table
go_df.to_csv(output_file_all, index=False)
print(f"[INFO] Full GO term table saved to: {output_file_all}")

# Create filtered table with highest GO_level per transcript
# For each transcript, keep rows where GO_level equals the max GO_level for that transcript
max_level_per_transcript = go_df.groupby("Transcript_ID")["GO_level"].transform("max")
filtered_go_df = go_df[go_df["GO_level"] == max_level_per_transcript].copy()

# Optional: sort nicely
filtered_go_df = filtered_go_df.sort_values(
    by=["Transcript_ID", "GO_level", "GO_term"],
    ascending=[True, False, True]
)

filtered_go_df.to_csv(output_file_filtered, index=False)
print(f"[INFO] Filtered GO table saved to: {output_file_filtered}")

# Summary
summary_all = (
    go_df.groupby("Ontology")["GO_term"]
    .nunique()
    .reset_index(name="GO_Term_Count")
)

summary_filtered = (
    filtered_go_df.groupby("Ontology")["GO_term"]
    .nunique()
    .reset_index(name="Filtered_GO_Term_Count")
)

print("\n[INFO] Unique GO term counts by ontology (full table):")
print(summary_all.to_string(index=False))

print("\n[INFO] Unique GO term counts by ontology (filtered table):")
print(summary_filtered.to_string(index=False))

print(f"\n[INFO] Transcripts without GO terms: {skipped_transcripts}")
print(f"[INFO] Unique transcripts with GO terms (full): {go_df['Transcript_ID'].nunique()}")
print(f"[INFO] Unique GO terms total (full): {go_df['GO_term'].nunique()}")
print(f"[INFO] Rows in filtered table: {len(filtered_go_df)}")
print(f"[INFO] Unique transcripts in filtered table: {filtered_go_df['Transcript_ID'].nunique()}")

# Output (for GFF3 annotation)
output_file = "GO_annotation.txt"

# Group by transcript and combine GO terms
gff_df = (
    filtered_go_df.groupby("Transcript_ID")["GO_term"]
    .apply(lambda terms: ",".join(sorted(set(terms))))
    .reset_index()
)

# Save as tab-separated file (no header if prefer)
gff_df.to_csv(output_file, sep="\t", index=False, header=False)

print(f"[INFO] GFF3 annotation table saved to: {output_file}")
