# 05_export_cytoscape_graphml.py
# Exports a .graphml file with all node attributes including 'origin'
# (broth_only / isolate_only / shared_broth_isolate)
# Load this in Cytoscape → filter on 'origin' column to remove broth nodes
# Then use Style panel for pie chart coloring with all 11 isolates

import pandas as pd
import networkx as nx
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = REPO_ROOT / "results" / "msms_analysis_output"

# Load your filtered (broth-subtracted) nodes — these already have the origin column
nodes = pd.read_csv(OUT_DIR / "nodes_isolate_only.tsv", sep="\t", low_memory=False)
edges = pd.read_csv(OUT_DIR / "edges_isolate_only.tsv", sep="\t")

# Also load isolate-only nodes to get that flag
isolate_only_ids = set(
    pd.read_csv(OUT_DIR / "nodes_isolate_only.tsv", sep="\t")["cluster index"]
)
nodes["isolate_only"] = nodes["cluster index"].isin(isolate_only_ids)

# Add a clean single-column annotation status for Cytoscape coloring.
# Flags are informational only (colors nodes for review) -- nothing is dropped.
# 300 ppm sits in this dataset's natural gap between real matches (<~50 ppm)
# and clearly-wrong ones (>~300 ppm), so plausible analogue/adduct hits aren't flagged.
PPM_FLAG_LIMIT = 300

has_match = nodes["Compound_Name"].notna() & (nodes["Compound_Name"] != "")
flag_ion_mobility = has_match & nodes.get("LibraryName", pd.Series("", index=nodes.index)).str.contains("ION-MOBILITY", na=False)
flag_high_ppm = has_match & (nodes["MZErrorPPM"].abs() > PPM_FLAG_LIMIT)
flagged = flag_ion_mobility | flag_high_ppm

nodes["flag_ion_mobility"] = flag_ion_mobility
nodes["flag_high_ppm"] = flag_high_ppm
nodes["annotation_status"] = "unannotated"
nodes.loc[has_match & ~flagged, "annotation_status"] = "valid_match"
nodes.loc[flagged, "annotation_status"] = "flagged"

print(f"Library hits flagged for review (kept in export, not removed):")
print(f"  Ion-mobility library source:   {flag_ion_mobility.sum()}")
print(f"  Mass error > {PPM_FLAG_LIMIT} ppm:          {flag_high_ppm.sum()}")
print(f"  Total flagged (either reason): {flagged.sum()}")

# Short compound label (for node labels in Cytoscape)
nodes["label"] = nodes["Compound_Name"].fillna("").apply(
    lambda x: x[:25] if isinstance(x, str) and x else ""
)

# Build graph
G = nx.Graph()

ISOLATE_COLS = [c for c in nodes.columns
                if c.startswith("ATTRIBUTE_ORGANISM:GNPSGROUP:")
                and c != "ATTRIBUTE_ORGANISM:GNPSGROUP:nan"]

# Add nodes with selected attributes (graphml can't handle list columns)
KEEP_ATTRS = [
    "cluster index", "number of spectra", "parent mass", "precursor charge",
    "RTMean", "component", "origin", "isolate_only", "annotation_status",
    "flag_ion_mobility", "flag_high_ppm",
    "Compound_Name", "label", "MQScore", "MZErrorPPM", "Adduct",
    "npclassifier_pathway", "npclassifier_class", "npclassifier_superclass",
    "Smiles", "InChIKey",
] + ISOLATE_COLS

KEEP_ATTRS = [c for c in KEEP_ATTRS if c in nodes.columns]

for _, row in nodes.iterrows():
    node_id = int(row["cluster index"])
    attrs = {}
    for col in KEEP_ATTRS:
        val = row[col]
        # graphml requires simple types — convert NaN to empty string
        if pd.isna(val):
            attrs[col.replace(":", "_").replace(" ", "_")] = ""
        elif isinstance(val, (int, float, bool, str)):
            attrs[col.replace(":", "_").replace(" ", "_")] = val
        else:
            attrs[col.replace(":", "_").replace(" ", "_")] = str(val)
    G.add_node(node_id, **attrs)

# Add edges
for _, row in edges.iterrows():
    n1, n2 = int(row["CLUSTERID1"]), int(row["CLUSTERID2"])
    if G.has_node(n1) and G.has_node(n2):
        G.add_edge(n1, n2,
                   cosine=float(row["Cosine"]),
                   delta_mz=float(row["DeltaMZ"]),
                   component=int(row["ComponentIndex"]) if "ComponentIndex" in row else -1)

# Export
out_path = OUT_DIR / "network_for_cytoscape.graphml"
nx.write_graphml(G, out_path)
print(f"Exported {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
print(f"Saved → {out_path}")
print()
print("In Cytoscape:")
print("  File → Import → Network from File → select network_for_cytoscape.graphml")
print("  File → Import → Table from File → (not needed, attributes are embedded)")
print()
print("To remove broth/shared nodes:")
print("  Select → Nodes → By Filter")
print("  Add filter: Column 'origin' IS 'isolate_only'  → Apply → Invert → Delete")
print("  OR: Column 'isolate_only' IS true → Apply → Invert → Delete")
print()
print("To highlight top components:")
print("  Select → Nodes → By Filter → Column 'component' IS <your_component_id>")
print("  Then: Edit → New Network from Selection")