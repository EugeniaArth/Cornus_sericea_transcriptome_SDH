from Bio import AlignIO
from collections import Counter
import math
import csv

# =========================
# USER SETTINGS
# =========================

alignment_file = "/Users/eugenianikonorova/Desktop/Work/Cornus_sericea_transcriptome_SDH/Files/Alignment/For_alignment_total_copy.fasta"


# Gallate-forming enzymes
positive_names = [
    "CsDQD/SDHc",
    "CsDQD/SDHd",
    "EcDQD/SDH3",
    "EcDQD/SDH2",
    "EgDQD/SDH3",
    "PgSDH3_1",
    "PgSDH3_2",
    "PgSDH3a_1",
    "PgSDH3a_2",
    "VvSDH4",
    "PgSDH4"
]

# Add true non-gallate enzymes here
# Strongly recommended: include Arabidopsis + other canonical SDHs
negative_names = [
    "AtSDH",
    "BrapaQDH",
    "BnapusQDH" ,
    "VvSDH2"
    # "Other_non_gallate_2",
]

arabidopsis_name = "AtSDH"  # exact/partial header match

# Known important A. thaliana residues from your figure
# Replace with actual Arabidopsis ungapped residue numbers if known
known_residues_substrate = [338, 381, 385, 422]
known_residues_cofactor = [483, 484, 485]
known_residues_all = sorted(set(known_residues_substrate + known_residues_cofactor))

# Hard filters
min_pos_conservation = 0.85
min_neg_conservation = 0.60
max_gap_fraction_pos = 0.20
max_gap_fraction_neg = 0.30
min_score_to_report = 7

# Neighborhood handling
cluster_window = 2   # collapse nearby hits +/- 2 columns

# =========================
# AA CHEMISTRY
# =========================

aa_class = {
    "A": "hydrophobic", "V": "hydrophobic", "I": "hydrophobic", "L": "hydrophobic",
    "M": "hydrophobic", "F": "aromatic", "W": "aromatic", "Y": "aromatic",
    "S": "polar", "T": "polar", "N": "polar", "Q": "polar", "C": "polar",
    "K": "positive", "R": "positive", "H": "positive",
    "D": "negative", "E": "negative",
    "G": "special", "P": "special",
    "-": "gap"
}

# BLOSUM-like rough severity categories without external dependencies
# higher = more radical
def substitution_severity(a, b):
    if a == b:
        return 0
    ca = aa_class.get(a, "other")
    cb = aa_class.get(b, "other")

    if a == "-" or b == "-":
        return 3
    if ca == cb:
        return 1
    if {"positive", "negative"} == {ca, cb}:
        return 4
    if "special" in {ca, cb}:
        return 3
    if {"hydrophobic", "aromatic"} == {ca, cb}:
        return 2
    return 3

def chemical_difference(a, b):
    return aa_class.get(a, "other") != aa_class.get(b, "other")

# =========================
# HELPERS
# =========================

def match_records(alignment, name_patterns):
    matched = []
    for rec in alignment:
        rid = rec.id
        if any(pat in rid for pat in name_patterns):
            matched.append(rec)
    return matched

def get_consensus_and_conservation(column):
    residues = [x for x in column if x != "-"]
    if not residues:
        return None, 0.0, Counter()
    cnt = Counter(residues)
    aa, n = cnt.most_common(1)[0]
    return aa, n / len(residues), cnt

def gap_fraction(column):
    if not column:
        return 1.0
    return sum(1 for x in column if x == "-") / len(column)

def shannon_entropy(column):
    residues = [x for x in column if x != "-"]
    if not residues:
        return 0.0
    cnt = Counter(residues)
    total = len(residues)
    ent = 0.0
    for n in cnt.values():
        p = n / total
        ent -= p * math.log2(p)
    return ent

def alignment_to_ungapped_map(seq):
    """
    alignment index -> ungapped residue number (1-based), or None if gap in this seq
    """
    mapping = {}
    pos = 0
    for i, aa in enumerate(str(seq)):
        if aa != "-":
            pos += 1
            mapping[i] = pos
        else:
            mapping[i] = None
    return mapping

def distance_to_known(ungapped_pos, known_positions):
    if ungapped_pos is None or not known_positions:
        return None
    return min(abs(ungapped_pos - k) for k in known_positions)

def neighborhood_bonus(dist):
    if dist is None:
        return 0
    if dist == 0:
        return 3
    if dist <= 2:
        return 3
    if dist <= 5:
        return 2
    if dist <= 8:
        return 1
    return 0

def is_indel_like(pos_col, neg_col):
    """
    Reward positions where positives and negatives differ in gap pattern.
    """
    pos_gap = gap_fraction(pos_col)
    neg_gap = gap_fraction(neg_col)
    return abs(pos_gap - neg_gap) >= 0.5

def collapse_hits(sorted_hits, window=2):
    """
    Keep only best-scoring hit in local neighborhoods.
    """
    kept = []
    for hit in sorted_hits:
        pos = hit["alignment_pos"]
        too_close = False
        for k in kept:
            if abs(pos - k["alignment_pos"]) <= window:
                too_close = True
                break
        if not too_close:
            kept.append(hit)
    return kept

# =========================
# LOAD ALIGNMENT
# =========================

alignment = AlignIO.read(alignment_file, "fasta")

positive_records = match_records(alignment, positive_names)
negative_records = match_records(alignment, negative_names)
arabidopsis_records = match_records(alignment, [arabidopsis_name])

if not positive_records:
    raise ValueError("No positive-group sequences matched.")
if not negative_records:
    raise ValueError("No negative-group sequences matched. Add non-gallate sequences.")
if len(arabidopsis_records) != 1:
    raise ValueError("Arabidopsis reference was not found uniquely.")

arabidopsis = arabidopsis_records[0]
at_map = alignment_to_ungapped_map(arabidopsis.seq)

aln_len = alignment.get_alignment_length()

# =========================
# SCAN COLUMNS
# =========================

hits = []

for i in range(aln_len):
    pos_col = [str(rec.seq[i]) for rec in positive_records]
    neg_col = [str(rec.seq[i]) for rec in negative_records]
    at_aa = str(arabidopsis.seq[i])

    pos_gap = gap_fraction(pos_col)
    neg_gap = gap_fraction(neg_col)

    # Hard filter 1: too many gaps
    if pos_gap > max_gap_fraction_pos:
        continue
    if neg_gap > max_gap_fraction_neg:
        continue

    pos_cons, pos_conv, pos_counts = get_consensus_and_conservation(pos_col)
    neg_cons, neg_conv, neg_counts = get_consensus_and_conservation(neg_col)

    if pos_cons is None or neg_cons is None:
        continue

    # Hard filter 2: positives must be strongly conserved
    if pos_conv < min_pos_conservation:
        continue

    # Hard filter 3: negatives should not be random noise
    if neg_conv < min_neg_conservation:
        continue

    # Hard filter 4: if pos and neg consensus same, not interesting
    if pos_cons == neg_cons:
        continue

    # Hard filter 5: skip columns where Arabidopsis itself has a gap
    if at_aa == "-":
        continue

    score = 0
    reasons = []

    # Strong positive-group conservation
    if pos_conv >= 0.95:
        score += 3
        reasons.append("pos_cons>=0.95")
    elif pos_conv >= 0.85:
        score += 2
        reasons.append("pos_cons>=0.85")

    # Strong negative-group conservation
    if neg_conv >= 0.90:
        score += 2
        reasons.append("neg_cons>=0.90")
    elif neg_conv >= 0.75:
        score += 1
        reasons.append("neg_cons>=0.75")

    # Positive consensus differs from Arabidopsis
    if pos_cons != at_aa:
        score += 2
        reasons.append("diff_from_At")

    # Positive consensus differs from negative consensus
    if pos_cons != neg_cons:
        score += 3
        reasons.append("diff_pos_vs_neg")

    # Chemistry
    if chemical_difference(pos_cons, at_aa):
        score += 1
        reasons.append("chem_diff_vs_At")

    sev_at = substitution_severity(pos_cons, at_aa)
    if sev_at >= 3:
        score += 2
        reasons.append("radical_vs_At")
    elif sev_at == 2:
        score += 1
        reasons.append("moderate_vs_At")

    sev_neg = substitution_severity(pos_cons, neg_cons)
    if sev_neg >= 3:
        score += 2
        reasons.append("radical_pos_vs_neg")
    elif sev_neg == 2:
        score += 1
        reasons.append("moderate_pos_vs_neg")

    # Downweight trivial residues
    if pos_cons in {"G", "A"}:
        score -= 1
        reasons.append("trivial_pos_residue")

    # Entropy bonus/penalty
    pos_ent = shannon_entropy(pos_col)
    neg_ent = shannon_entropy(neg_col)

    if pos_ent < 0.3:
        score += 1
        reasons.append("low_pos_entropy")
    if neg_ent > 1.0:
        score -= 1
        reasons.append("high_neg_entropy")

    # Indel-like difference
    if is_indel_like(pos_col, neg_col):
        score += 1
        reasons.append("indel_pattern")

    # Proximity to known Arabidopsis residues
    at_pos = at_map[i]
    dist = distance_to_known(at_pos, known_residues_all)
    nb = neighborhood_bonus(dist)
    if nb:
        score += nb
        reasons.append(f"near_known_{dist}")

    # Optional: require at least one of the biologically meaningful criteria
    meaningful = (
        (pos_cons != neg_cons) and
        (pos_cons != at_aa or dist is not None and dist <= 8)
    )
    if not meaningful:
        continue

    if score >= min_score_to_report:
        hits.append({
            "alignment_pos": i + 1,
            "arabidopsis_pos": at_pos,
            "at_residue": at_aa,
            "positive_consensus": pos_cons,
            "positive_conservation": round(pos_conv, 3),
            "negative_consensus": neg_cons,
            "negative_conservation": round(neg_conv, 3),
            "score": score,
            "distance_to_known": dist,
            "reasons": ",".join(reasons),
            "pos_column": "".join(pos_col),
            "neg_column": "".join(neg_col),
        })

# =========================
# SORT AND COLLAPSE
# =========================

hits.sort(
    key=lambda x: (
        -x["score"],
        9999 if x["distance_to_known"] is None else x["distance_to_known"],
        x["alignment_pos"]
    )
)

collapsed_hits = collapse_hits(hits, window=cluster_window)

# =========================
# OUTPUT
# =========================

print(f"\nRaw hits: {len(hits)}")
print(f"Collapsed hits: {len(collapsed_hits)}\n")

header = [
    "alignment_pos", "arabidopsis_pos", "at_residue",
    "positive_consensus", "positive_conservation",
    "negative_consensus", "negative_conservation",
    "score", "distance_to_known", "reasons"
]

print("\t".join(header))
for h in collapsed_hits:
    print("\t".join(str(h[k]) for k in header))

with open("candidate_positions.tsv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "alignment_pos", "arabidopsis_pos", "at_residue",
        "positive_consensus", "positive_conservation",
        "negative_consensus", "negative_conservation",
        "score", "distance_to_known", "reasons",
        "pos_column", "neg_column"
    ], delimiter="\t")
    writer.writeheader()
    writer.writerows(collapsed_hits)

print("\nSaved: candidate_positions.tsv")
