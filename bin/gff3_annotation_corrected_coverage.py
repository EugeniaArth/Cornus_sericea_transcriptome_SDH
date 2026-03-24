import pandas as pd
import re

# -----------------------------
# ID normalization
# removes the variable coverage part:
# NODE_10_length_27663_cov_217.671471_g1_i3 -> NODE_10_length_27663_g1_i3
# -----------------------------
def normalize_id(x):
    x = str(x).strip()
    return re.sub(r'_cov_[^_]+', '', x)

# -----------------------------
# Load RefSeq annotations (highest priority)
# -----------------------------
refseq = pd.read_csv(
    'Files/refseq_best_hits.txt',
    sep='\t',
    header=None,
    names=['ID', 'Hit', 'E-value', 'Score', 'Description']
)
refseq['norm_ID'] = refseq['ID'].apply(normalize_id)
refseq_dict = refseq.drop_duplicates('norm_ID').set_index('norm_ID').to_dict(orient='index')

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
    description = str(description)

    gene_match = re.search(r'GN=([^ ]+)', description)
    gene_name = gene_match.group(1) if gene_match else 'NA'

    product_match = re.search(r'sp\|[^|]+\|[^ ]+\s+(.+?)\s+OS=', description)
    product = product_match.group(1) if product_match else 'NA'

    uniprot_match = re.match(r'(sp\|[^ ]+)', description)
    uniprot_hit = uniprot_match.group(1) if uniprot_match else 'NA'

    return gene_name, product, uniprot_hit

uniprot['norm_ID'] = uniprot['ID'].apply(normalize_id)

uniprot_dict = {}
for _, row in uniprot.drop_duplicates('norm_ID').iterrows():
    gene_name, product, uniprot_hit = parse_uniprot(row['Description'])
    uniprot_dict[row['norm_ID']] = {
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
nr['norm_ID'] = nr['ID'].apply(normalize_id)
nr_dict = nr.drop_duplicates('norm_ID').set_index('norm_ID').to_dict(orient='index')

# -----------------------------
# Load GO annotations (only from filtered GO file)
# -----------------------------
go_annots = pd.read_csv(
    'GO_annotation.txt',
    sep='\t',
    header=None,
    names=['ID', 'GO_terms']
)
go_annots['norm_ID'] = go_annots['ID'].apply(normalize_id)
go_dict = go_annots.drop_duplicates('norm_ID').set_index('norm_ID')['GO_terms'].to_dict()

# -----------------------------
# Load KEGG annotations per transcript
# -----------------------------
kegg_df = pd.read_csv('Files/KEGG_annotation_per_transcript.csv', sep=None, engine='python')
kegg_df.columns = kegg_df.columns.str.strip()
print("KEGG columns:", list(kegg_df.columns))

if 'Transcript_ID' not in kegg_df.columns:
    raise ValueError(f"Expected 'Transcript_ID' column, found: {kegg_df.columns.tolist()}")

for col in ['KO_IDs', 'Pathway_IDs', 'Module_IDs', 'Reaction_IDs']:
    if col not in kegg_df.columns:
        kegg_df[col] = ''

kegg_df = kegg_df.fillna('')
kegg_df['norm_ID'] = kegg_df['Transcript_ID'].apply(normalize_id)

def clean_kegg_ko(value):
    if not str(value).strip():
        return ''
    return ','.join(
        item.strip().replace('ko:', '')
        for item in str(value).split(',')
        if item.strip()
    )

kegg_df['KO_IDs'] = kegg_df['KO_IDs'].apply(clean_kegg_ko)
kegg_dict = kegg_df.drop_duplicates('norm_ID').set_index('norm_ID').to_dict(orient='index')

# -----------------------------
# Minimal checks
# -----------------------------
print("RefSeq entries:", len(refseq_dict))
print("UniProt entries:", len(uniprot_dict))
print("nr entries:", len(nr_dict))

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
            raw_lookup_id = attr_fields.get('ID', '').strip()
        else:
            raw_lookup_id = attr_fields.get('Parent', '').strip()

        lookup_id = normalize_id(raw_lookup_id)

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
                unannotated_ids.append(raw_lookup_id)

            go_term = go_dict.get(lookup_id, '')
            if pd.notna(go_term) and str(go_term).strip():
                go_count += 1

            if lookup_id in kegg_dict and any(
                str(kegg_dict[lookup_id].get(col, '')).strip()
                for col in ['KO_IDs', 'Pathway_IDs', 'Module_IDs', 'Reaction_IDs']
            ):
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
        # GO annotation
        # -----------------------------
        go_term = go_dict.get(lookup_id, '')
        if pd.notna(go_term) and str(go_term).strip():
            attr_fields['Ontology_term'] = str(go_term)

        # -----------------------------
        # KEGG annotation
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

        parts[8] = ';'.join(f'{k}={v}' for k, v in attr_fields.items())
        outfile.write('\t'.join(parts) + '\n')

# -----------------------------
# Debug overlap after normalization
# -----------------------------
norm_nr_keys = set(nr_dict.keys())
print("Normalized overlap between GFF mRNA and nr:", len(seen_transcripts & norm_nr_keys))

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

with open("unannotated_transcripts.txt", "w") as f:
    for uid in unannotated_ids:
        f.write(uid + "\n")

print("\nUnannotated transcript IDs saved to: unannotated_transcripts.txt")
