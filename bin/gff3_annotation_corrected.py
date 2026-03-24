import pandas as pd
import re

# -----------------------------
# Load RefSeq annotations (highest priority)
# -----------------------------
refseq = pd.read_csv(
    'Files/refseq_best_hits.txt',
    sep='\t',
    header=None,
    names=['ID', 'Hit', 'E-value', 'Score', 'Description']
)
refseq_dict = refseq.set_index('ID').to_dict(orient='index')

# -----------------------------
# Load UniProt annotations (second priority)
# -----------------------------
uniprot = pd.read_csv(
    'Files/uniprot_best_hits.txt',
    sep='\t',
    header=None,
    names=['ID', 'Hit', 'E-value', 'Score', 'Description']
)

def parse_uniprot(description):
    gene_match = re.search(r'GN=([^ ]+)', str(description))
    gene_name = gene_match.group(1) if gene_match else 'NA'

    product_match = re.search(r'sp\|[^|]+\|[^ ]+\s+(.+?)\s+OS=', str(description))
    product = product_match.group(1) if product_match else 'NA'

    uniprot_match = re.match(r'(sp\|[^ ]+)', str(description))
    uniprot_hit = uniprot_match.group(1) if uniprot_match else 'NA'

    return gene_name, product, uniprot_hit

uniprot_dict = {}
for _, row in uniprot.iterrows():
    gene_name, product, uniprot_hit = parse_uniprot(row['Description'])
    uniprot_dict[row['ID']] = {
        'gene_name': gene_name,
        'product': product,
        'uniprot_hit': uniprot_hit
    }

# -----------------------------
# Load NCBI nr annotations (third priority)
# -----------------------------
nr = pd.read_csv(
    'Files/nr_best_hits.txt',
    sep='\t',
    header=None,
    names=['ID', 'Hit', 'E-value', 'Score', 'Description']
)
nr_dict = nr.set_index('ID').to_dict(orient='index')

# -----------------------------
# Load GO annotations (only from filtered GO file)
# -----------------------------
go_annots = pd.read_csv(
    'Files/GO_annotation.txt',
    sep='\t',
    header=None,
    names=['ID', 'GO_terms']
)
go_dict = go_annots.set_index('ID')['GO_terms'].to_dict()

# -----------------------------
# Load KEGG annotations per transcript
# -----------------------------
kegg_df = pd.read_csv('KEGG_annotation_per_transcript.csv', sep=None, engine='python')

# Clean column names (remove spaces/BOM issues)
kegg_df.columns = kegg_df.columns.str.strip()

# Debug: print detected columns
print("KEGG columns:", list(kegg_df.columns))

# Ensure correct column name
if 'Transcript_ID' not in kegg_df.columns:
    raise ValueError(f"Expected 'Transcript_ID' column, found: {kegg_df.columns.tolist()}")

for col in ['KO_IDs', 'Pathway_IDs', 'Module_IDs', 'Reaction_IDs']:
    if col not in kegg_df.columns:
        kegg_df[col] = ''

kegg_df = kegg_df.fillna('')

def clean_kegg_ko(value):
    if not str(value).strip():
        return ''
    return ','.join(
        item.strip().replace('ko:', '')
        for item in str(value).split(',')
        if item.strip()
    )

kegg_df['KO_IDs'] = kegg_df['KO_IDs'].apply(clean_kegg_ko)

kegg_dict = kegg_df.set_index('Transcript_ID').to_dict(orient='index')

# -----------------------------
# Counters
# -----------------------------
refseq_count = 0
uniprot_count = 0
nr_count = 0
unannotated_count = 0
go_count = 0
kegg_count = 0
seen_transcripts = set()
unannotated_ids = []

# -----------------------------
# Process GFF3
# -----------------------------
with open('Files/final.clust_annotation_longest_iso.gff3', 'r') as infile, \
     open('annotated_final.gff3', 'w') as outfile:

    for line in infile:
        if line.startswith('#') or line.strip() == '':
            outfile.write(line)
            continue

        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9:
            outfile.write(line)
            continue

        feature_type = parts[2]
        attributes = parts[8]

        attr_fields = {}
        for field in attributes.split(';'):
            if '=' in field:
                key, value = field.split('=', 1)
                attr_fields[key] = value

        # Use ID for gene/mRNA, Parent for child features
        if feature_type in ('gene', 'mRNA'):
            lookup_id = attr_fields.get('ID', '').strip()
        else:
            lookup_id = attr_fields.get('Parent', '').strip()

        # Count only once per transcript (mRNA features only)
        if feature_type == 'mRNA' and lookup_id not in seen_transcripts:
            seen_transcripts.add(lookup_id)

            if lookup_id in refseq_dict:
                refseq_count += 1
            elif lookup_id in uniprot_dict:
                uniprot_count += 1
            elif lookup_id in nr_dict:
                nr_count += 1
            else:
                unannotated_count += 1
                unannotated_ids.append(lookup_id)

            go_term = go_dict.get(lookup_id, '')
            if pd.notna(go_term) and str(go_term).strip():
                go_count += 1

            if lookup_id in kegg_dict and any(str(kegg_dict[lookup_id].get(col, '')).strip() for col in ['KO_IDs', 'Pathway_IDs', 'Module_IDs', 'Reaction_IDs']):
                kegg_count += 1

        # -----------------------------
        # Functional annotation priority:
        # RefSeq > UniProt > nr
        # -----------------------------
        if lookup_id in refseq_dict:
            attr_fields['Name'] = refseq_dict[lookup_id]['Hit']
            attr_fields['product'] = refseq_dict[lookup_id]['Description']
            attr_fields['RefSeq_best_hit'] = refseq_dict[lookup_id]['Hit']

        elif lookup_id in uniprot_dict:
            attr_fields['Name'] = uniprot_dict[lookup_id]['gene_name']
            attr_fields['product'] = uniprot_dict[lookup_id]['product']
            attr_fields['UniProt_best_hit'] = uniprot_dict[lookup_id]['uniprot_hit']

        elif lookup_id in nr_dict:
            attr_fields['Name'] = nr_dict[lookup_id]['Hit']
            attr_fields['product'] = nr_dict[lookup_id]['Description']
            attr_fields['NR_best_hit'] = nr_dict[lookup_id]['Hit']

        # -----------------------------
        # GO annotation (only from GO_annotation.txt)
        # -----------------------------
        go_term = go_dict.get(lookup_id, '')
        if pd.notna(go_term) and str(go_term).strip():
            attr_fields['Ontology_term'] = str(go_term)

        # -----------------------------
        # KEGG annotation (from KEGG_annotation_per_transcript.csv)
        # -----------------------------
        if lookup_id in kegg_dict:
            kegg_row = kegg_dict[lookup_id]

            if str(kegg_row.get('KO_IDs', '')).strip():
                attr_fields['KEGG_KO'] = str(kegg_row['KO_IDs'])

            if str(kegg_row.get('Pathway_IDs', '')).strip():
                attr_fields['KEGG_Pathway'] = str(kegg_row['Pathway_IDs'])

            if str(kegg_row.get('Module_IDs', '')).strip():
                attr_fields['KEGG_Module'] = str(kegg_row['Module_IDs'])

            if str(kegg_row.get('Reaction_IDs', '')).strip():
                attr_fields['KEGG_Reaction'] = str(kegg_row['Reaction_IDs'])

        # Rebuild attributes
        parts[8] = ';'.join(f'{k}={v}' for k, v in attr_fields.items())
        outfile.write('\t'.join(parts) + '\n')

# -----------------------------
# Print summary
# -----------------------------
print("Annotated GFF3 file created: annotated_final.gff3")
print("\nAnnotation Summary (mRNA features only):")
print(f"RefSeq annotations: {refseq_count}")
print(f"UniProt annotations: {uniprot_count}")
print(f"NCBI nr annotations: {nr_count}")
print(f"GO annotations: {go_count}")
print(f"KEGG annotations: {kegg_count}")
print(f"Unannotated transcripts: {unannotated_count}")

# Save unannotated IDs
with open("unannotated_transcripts.txt", "w") as f:
    for uid in unannotated_ids:
        f.write(uid + "\n")

print("\nUnannotated transcript IDs saved to: unannotated_transcripts.txt")
