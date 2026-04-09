from Bio import AlignIO
from collections import Counter

# =========================
# USER SETTINGS
# =========================

alignment_file = "Files/Alignment/claster_clustalO.fa"

# names must match FASTA headers (partial match works)
gallate_group_names = [
    "CsDQD", "EcDQD", "EgDQD",
    "PgSDH3_1", "PgSDH3_2",
    "PgSDH3a_1", "PgSDH3a_2",
    "VvSDH4", "PgSDH4"
]

reference_name = "AtSDH"  # adjust if needed

conservation_threshold = 0.8  # 80% conserved

# =========================
# AA GROUPS FOR CHEMISTRY
# =========================

aa_groups = {
    "nonpolar": set("AVLIMFWP"),
    "polar": set("STNQYC"),
    "positive": set("KRH"),
    "negative": set("DE"),
    "special": set("G")
}

def get_group(aa):
    for group, members in aa_groups.items():
        if aa in members:
            return group
    return "other"

# =========================
# LOAD ALIGNMENT
# =========================

alignment = AlignIO.read(alignment_file, "fasta")

# separate sequences
gallate_seqs = []
reference_seq = None

for record in alignment:
    name = record.id

    if any(x in name for x in gallate_group_names):
        gallate_seqs.append(record.seq)

    if reference_name in name:
        reference_seq = record.seq

if reference_seq is None:
    raise ValueError("Reference sequence not found!")

# =========================
# SCORING FUNCTION
# =========================

def score_position(column, ref_aa):
    score = 0

    # remove gaps
    residues = [aa for aa in column if aa != "-"]

    if len(residues) == 0:
        return None

    counts = Counter(residues)
    most_common_aa, freq = counts.most_common(1)[0]
    conservation = freq / len(residues)

    # +2 conserved in gallate group
    if conservation >= conservation_threshold:
        score += 2

    # +2 different from reference
    if most_common_aa != ref_aa:
        score += 2

    # +1 chemical difference
    if get_group(most_common_aa) != get_group(ref_aa):
        score += 1

    # +1 not trivial residue
    if most_common_aa not in ["G", "A"]:
        score += 1

    return {
        "score": score,
        "consensus": most_common_aa,
        "conservation": round(conservation, 2)
    }

# =========================
# MAIN LOOP
# =========================

results = []

alignment_length = alignment.get_alignment_length()

for i in range(alignment_length):
    column = [seq[i] for seq in gallate_seqs]
    ref_aa = reference_seq[i]

    if ref_aa == "-":
        continue

    result = score_position(column, ref_aa)

    if result is None:
        continue

    if result["score"] >= 4:  # threshold for reporting
        results.append({
            "position": i + 1,
            "ref": ref_aa,
            "consensus": result["consensus"],
            "score": result["score"],
            "conservation": result["conservation"]
        })

# =========================
# OUTPUT
# =========================

results = sorted(results, key=lambda x: x["score"], reverse=True)

print("\nTop candidate positions:\n")
print("Pos\tRef\tGallate\tScore\tConservation")

for r in results:
    print(f"{r['position']}\t{r['ref']}\t{r['consensus']}\t{r['score']}\t{r['conservation']}")
