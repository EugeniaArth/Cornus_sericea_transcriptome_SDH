from pathlib import Path
import pandas as pd


INPUT_FILE = "Files/eggnog_annotations.tsv"
INTERMEDIATE_FILE = "Files/KEGG/KEGG_annotation_per_transcript.csv"
OUTPUT_FILE = "Files/KEGG_annotation_description.csv"

KO_FILE = "Files/KEGG/ko_list.txt"
PATHWAY_FILE = "Files/KEGG/pathway_list.txt"
REACTION_FILE = "Files/KEGG/reaction_list.txt"
MODULE_FILE = "Files/KEGG/kegg_modules.txt"


# Helpers

def load_lookup(filename):
    """Load 2-column tab/space-delimited lookup file into a dict."""
    lookup = {}
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 2:
                parts = line.split()

            if len(parts) >= 2:
                lookup[parts[0].strip()] = " ".join(parts[1:]).strip()
    return lookup


def load_module_dict(filename):
    """Parse KEGG module file and return module_id -> description/group/prefix."""
    modules = {}
    current_prefix = ""
    current_group = ""

    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("B") and "<b>" in line and "</b>" in line:
                current_prefix = line.split("<b>")[-1].split("</b>")[0].strip()

            elif line.startswith("C"):
                current_group = line[1:].strip()

            elif line.startswith("D"):
                parts = line.split()
                if len(parts) >= 3:
                    module_id = parts[1]
                    rest = " ".join(parts[2:])
                    description = rest.split("[")[0].strip()
                    modules[module_id] = {
                        "description": description,
                        "group": current_group,
                        "prefix": current_prefix,
                    }
    return modules


def split_ids(value):
    """Split comma-separated annotation field safely."""
    if pd.isna(value) or not str(value).strip():
        return []
    return [x.strip() for x in str(value).split(",") if x.strip()]


def map_ids(ids, lookup, transform=None, keep_if=None):
    """Map IDs through a lookup dict with optional transform/filter."""
    result = []
    for item in ids:
        key = transform(item) if transform else item
        if keep_if and not keep_if(item):
            continue
        result.append(lookup.get(key, key))
    return result


def map_modules(ids, module_dict):
    """Return descriptions, groups, and prefixes for module IDs."""
    desc, groups, prefixes = [], [], []

    for m in ids:
        info = module_dict.get(m)
        if info:
            desc.append(info["description"])
            groups.append(info["group"])
            prefixes.append(info["prefix"])
        else:
            desc.append(m)
            groups.append("")
            prefixes.append("")

    return pd.Series({
        "Module_Descriptions": ";".join(desc),
        "Module_Groups": ";".join(groups),
        "Module_Prefixes": ";".join(prefixes),
    })


# Extract KEGG fields

df = pd.read_csv(INPUT_FILE, sep="\t", header=None, low_memory=False)

df_kegg = df[[0, 11, 12, 13, 14]].copy()
df_kegg.columns = [
    "Transcript_ID",
    "KO_IDs",
    "Pathway_IDs",
    "Module_IDs",
    "Reaction_IDs",
]

df_kegg = df_kegg.replace("-", pd.NA)
df_kegg = df_kegg.dropna(
    how="all",
    subset=["KO_IDs", "Pathway_IDs", "Module_IDs", "Reaction_IDs"]
)

df_kegg.to_csv(INTERMEDIATE_FILE, index=False)
print(f"[✓] KEGG annotations saved to: {INTERMEDIATE_FILE}")
print(df_kegg.head())


# Load lookups

ko_dict = load_lookup(KO_FILE)
pathway_dict = load_lookup(PATHWAY_FILE)
reaction_dict = load_lookup(REACTION_FILE)
module_dict = load_module_dict(MODULE_FILE)


# Add descriptions

df_ann = pd.read_csv(INTERMEDIATE_FILE)

df_ann["KO_IDs"] = df_ann["KO_IDs"].apply(
    lambda x: ";".join(
        map_ids(
            split_ids(x),
            ko_dict,
            transform=lambda s: s.replace("ko:", "")
        )
    )
)

df_ann["Pathway_IDs"] = df_ann["Pathway_IDs"].apply(
    lambda x: ";".join(
        map_ids(
            split_ids(x),
            pathway_dict,
            keep_if=lambda s: s.startswith("map")
        )
    )
)

df_ann["Reaction_IDs"] = df_ann["Reaction_IDs"].apply(
    lambda x: ";".join(map_ids(split_ids(x), reaction_dict))
)

module_info = df_ann["Module_IDs"].apply(lambda x: map_modules(split_ids(x), module_dict))
df_out = pd.concat([df_ann[["Transcript_ID", "KO_IDs", "Pathway_IDs", "Reaction_IDs"]], module_info], axis=1)

df_out.to_csv(OUTPUT_FILE, index=False)


# Summary stats

total = len(df_out)
with_module = (df_out["Module_Descriptions"] != "").sum()
with_pathway = (df_out["Pathway_IDs"] != "").sum()

print(f"\n {OUTPUT_FILE} created.")
print("📊 Annotation Summary:")
print(f"  Total Transcript_IDs         : {total}")
print(f"  With KEGG Modules            : {with_module} ({with_module / total:.1%})")
print(f"  Missing KEGG Modules         : {total - with_module}")
print(f"  With KEGG Pathways           : {with_pathway} ({with_pathway / total:.1%})")
print(f"  Missing KEGG Pathways        : {total - with_pathway}")
