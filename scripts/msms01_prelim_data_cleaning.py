# =============================================================================
# 00_gnps_cmn_analysis.py
# MicroChemEco Lab | R. Wang
#
# Primary input: nf_output/networking/clustersummary_with_network.tsv
# file contains nodes + per-isolate counts + library hits + classes
# =============================================================================

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")





# ── 0. PATHS ──────────────────────────────────────────────────────────────────
REPO_ROOT    = Path(__file__).resolve().parents[1]
NETWORK_DIR  = REPO_ROOT / "data" / "MCE-GNPS" / "nf_output" / "networking"
OUT_DIR      = REPO_ROOT / "results" / "msms_analysis_output"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Primary working file — has everything in one place
MAIN_FILE = NETWORK_DIR / "clustersummary_with_network.tsv"
EDGES_FILE = NETWORK_DIR / "filtered_pairs.tsv"   # or pairs_with_components.tsv





# ── 1. LOAD ───────────────────────────────────────────────────────────────────
print("Loading main node+network table...")
df = pd.read_csv(MAIN_FILE, sep="\t", low_memory=False)
print(f"  {len(df)} nodes, {len(df.columns)} columns")

edges = pd.read_csv(EDGES_FILE, sep="\t")
print(f"  {len(edges)} edges loaded")





# ── 2. DEFINE COLUMN GROUPS ───────────────────────────────────────────────────

BROTH_COL = "ATTRIBUTE_PERFILEGROUPING:GNPSGROUP:bkgrnd-marine-broth.mzML"

# 11 isolate columns (ORGANISM-level — aggregated per isolate)
ISOLATE_COLS = [
    "ATTRIBUTE_ORGANISM:GNPSGROUP:Vibrio_1A01",
    "ATTRIBUTE_ORGANISM:GNPSGROUP:Vibrio_6D03",
    "ATTRIBUTE_ORGANISM:GNPSGROUP:Phaeobacter_B2R04",
    "ATTRIBUTE_ORGANISM:GNPSGROUP:Shewanella_B2R08",
    "ATTRIBUTE_ORGANISM:GNPSGROUP:Roseovarius_B2R09",
    "ATTRIBUTE_ORGANISM:GNPSGROUP:Colwellia_C2M19",
    "ATTRIBUTE_ORGANISM:GNPSGROUP:Paracoccus_C2R07",
    "ATTRIBUTE_ORGANISM:GNPSGROUP:Salipiger_C3M06",
    "ATTRIBUTE_ORGANISM:GNPSGROUP:Acinetobacter_D3M06",
    "ATTRIBUTE_ORGANISM:GNPSGROUP:Arenibacter_D3M17",
    "ATTRIBUTE_ORGANISM:GNPSGROUP:Motilimonas_G1M02",
]

# Short display labels for plots (same order as above)
ISOLATE_LABELS = [
    "Vibrio 1A01", "Vibrio 6D03", "Phaeobacter B2R04",
    "Shewanella B2R08", "Roseovarius B2R09", "Colwellia C2M19",
    "Paracoccus C2R07", "Salipiger C3M06", "Acinetobacter D3M06",
    "Arenibacter D3M17", "Motilimonas G1M02",
]

# 11-color palette (colorblind-distinguishable)
COLORS = [
    "#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99",
    "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99"
]





# =============================================================================
# 3. SANITY CHECKS
# =============================================================================
print("\n" + "="*60)
print("SANITY CHECKS")
print("="*60)

total = len(df)
print(f"\nTotal nodes:   {total}")
print(f"Total edges:   {len(edges)}")

# Spectrum count distribution
print(f"\nSpectra per cluster:")
print(df["number of spectra"].describe().round(1).to_string())

# Charge distribution
print(f"\nCharge distribution:")
print(df["precursor charge"].value_counts().sort_index().to_string())
n_high_charge = (df["precursor charge"] >= 3).sum()
if n_high_charge:
    print(f"  ⚠ {n_high_charge} nodes with charge ≥3 — flag for review")

# Library match stats
has_match = df["Compound_Name"].notna() & (df["Compound_Name"] != "")
print(f"\nLibrary-matched nodes: {has_match.sum()} ({100*has_match.sum()/total:.1f}%)")
print(f"Unannotated nodes:     {(~has_match).sum()} ({100*(~has_match).sum()/total:.1f}%)")

if has_match.any():
    print(f"\nCosine score stats (library hits only):")
    print(df.loc[has_match, "MQScore"].describe().round(3).to_string())



# ── FLAG (do not remove): ion-mobility source or extreme mass error ─────────
# 300 ppm sits in the natural gap in this dataset's error distribution
# (real matches cluster under ~50 ppm, clearly-wrong ones start above ~300 ppm),
# so it catches obvious garbage without flagging plausible analogue/adduct hits.
PPM_FLAG_LIMIT = 300

if "MZErrorPPM" in df.columns:
    flag_ion_mobility = has_match & df.get("LibraryName", pd.Series("")).str.contains("ION-MOBILITY", na=False)
    flag_high_ppm = has_match & (df["MZErrorPPM"].abs() > PPM_FLAG_LIMIT)
    flagged = flag_ion_mobility | flag_high_ppm

    df["flag_ion_mobility"] = flag_ion_mobility
    df["flag_high_ppm"] = flag_high_ppm

    print(f"\n⚠ Library hits flagged for review (kept in data, not removed):")
    print(f"  Ion-mobility library source:   {flag_ion_mobility.sum()}")
    print(f"  Mass error > {PPM_FLAG_LIMIT} ppm:          {flag_high_ppm.sum()}")
    print(f"  Total flagged (either reason): {flagged.sum()}")

    df.loc[flagged, "match_quality"] = "FLAGGED"
    df.loc[has_match & ~flagged, "match_quality"] = "VALID"
    df.loc[~has_match, "match_quality"] = "UNANNOTATED"





# =============================================================================
# 4. BROTH SUBTRACTION (strict: drop anything present in broth at all, even
#    if it's also present in an isolate)
# =============================================================================
print("\n" + "="*60)
print("BROTH SUBTRACTION")
print("="*60)

broth_present   = df[BROTH_COL] > 0
isolate_present = df[ISOLATE_COLS].sum(axis=1) > 0

broth_only   = broth_present & ~isolate_present
isolate_only = ~broth_present & isolate_present
shared       = broth_present & isolate_present

df["origin"] = "unknown"
df.loc[broth_only,   "origin"] = "broth_only"
df.loc[isolate_only, "origin"] = "isolate_only"
df.loc[shared,       "origin"] = "shared_broth_isolate"

print(f"  Broth-only (REMOVED):            {broth_only.sum()}")
print(f"  Shared broth+isolate (REMOVED):  {shared.sum()}")
print(f"  Isolate-only (kept):             {isolate_only.sum()}")

df_f = df[isolate_only].copy()
print(f"\n  After broth subtraction: {len(df_f)} nodes retained ({len(df)-len(df_f)} removed)")

# Save isolate-only node table
df_f.to_csv(OUT_DIR / "nodes_isolate_only.tsv", sep="\t", index=False)
print(f"  Saved → analysis_output/nodes_isolate_only.tsv")

# Filter edges
kept_ids = set(df_f["cluster index"])
ef = edges[edges["CLUSTERID1"].isin(kept_ids) & edges["CLUSTERID2"].isin(kept_ids)].copy()
ef.to_csv(OUT_DIR / "edges_isolate_only.tsv", sep="\t", index=False)
print(f"  Saved → analysis_output/edges_isolate_only.tsv ({len(ef)} edges retained)")





# =============================================================================
# 5. LIBRARY HITS — CLEAN VALIDATION TABLE
# =============================================================================
print("\n" + "="*60)
print("LIBRARY HITS (all matches; see match_quality / flag_* columns)")
print("="*60)

valid_hits = df_f[
    df_f["Compound_Name"].notna() &
    (df_f["Compound_Name"] != "")
].copy()

n_flagged_hits = (valid_hits.get("match_quality", pd.Series()) == "FLAGGED").sum()
print(f"  Library hits on filtered nodes: {len(valid_hits)} ({n_flagged_hits} flagged)")

display_cols = [c for c in [
    "cluster index", "Compound_Name", "MQScore", "MZErrorPPM",
    "SharedPeaks", "Adduct", "npclassifier_pathway", "npclassifier_class", "origin",
    "match_quality", "flag_ion_mobility", "flag_high_ppm"
] if c in valid_hits.columns]

print(valid_hits[display_cols].sort_values("MQScore", ascending=False).to_string(index=False))
valid_hits.to_csv(OUT_DIR / "library_hits_clean.tsv", sep="\t", index=False)
print(f"\n  Saved → analysis_output/library_hits_clean.tsv")





# =============================================================================
# 6. CHEMICAL CLASS SUMMARY (NPClassifier)
# =============================================================================
print("\n" + "="*60)
print("NPClassifier CHEMICAL CLASS SUMMARY")
print("="*60)

if "npclassifier_pathway" in df_f.columns:
    pathway_counts = df_f["npclassifier_pathway"].dropna().value_counts()
    print("\nTop chemical pathways (valid annotated nodes):")
    print(pathway_counts.head(15).to_string())

    fig, ax = plt.subplots(figsize=(8, 5))
    pathway_counts.head(10).sort_values().plot(kind="barh", ax=ax, color="#1f78b4")
    ax.set_xlabel("Number of nodes")
    ax.set_title("Chemical class distribution (NPClassifier pathway)\nMarine isolates — broth-subtracted")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "chemical_class_barplot.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved → analysis_output/chemical_class_barplot.png")





# =============================================================================
# 7. ANNOTATION STATUS PIE CHART
# =============================================================================
flagged_n  = (df_f.get("match_quality", pd.Series()) == "FLAGGED").sum()
valid_n    = len(valid_hits) - flagged_n
unanno_n   = len(df_f) - valid_n - flagged_n

fig, ax = plt.subplots(figsize=(5, 5))
ax.pie(
    [valid_n, flagged_n, unanno_n],
    labels=[f"Library-matched\n(valid, n={valid_n})",
            f"Flagged hits\n(ion-mobility / high ppm, n={flagged_n})",
            f"Unannotated\n(n={unanno_n})"],
    colors=["#2196F3", "#FF9800", "#BDBDBD"],
    autopct="%1.1f%%", startangle=90,
    wedgeprops={"linewidth": 1, "edgecolor": "white"}
)
ax.set_title("Spectral annotation status\n(broth-subtracted nodes)", fontsize=11)
plt.tight_layout()
plt.savefig(OUT_DIR / "annotation_status_pie.png", dpi=150, bbox_inches="tight")
plt.close()
print("\nSaved → analysis_output/annotation_status_pie.png")





# =============================================================================
# 8. PER-ISOLATE FEATURE COUNT BARPLOT
# =============================================================================
isolate_totals = {}
for col, label in zip(ISOLATE_COLS, ISOLATE_LABELS):
    if col in df_f.columns:
        isolate_totals[label] = (df_f[col] > 0).sum()

fig, ax = plt.subplots(figsize=(10, 5))
bars = ax.bar(isolate_totals.keys(), isolate_totals.values(),
              color=COLORS[:len(isolate_totals)], edgecolor="white", linewidth=0.5)
ax.set_ylabel("Number of spectral features detected")
ax.set_title("Spectral features per isolate (broth-subtracted)", fontsize=12)
plt.xticks(rotation=35, ha="right", fontsize=9)
for bar, val in zip(bars, isolate_totals.values()):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
            str(val), ha="center", va="bottom", fontsize=8)
plt.tight_layout()
plt.savefig(OUT_DIR / "features_per_isolate.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved → analysis_output/features_per_isolate.png")





# =============================================================================
# 9. PIE CHART NETWORK (python-rendered, all 11 isolates)
# =============================================================================
# print("\n" + "="*60)
# print("BUILDING PIE CHART NETWORK")
# print("="*60)
#
# G = nx.Graph()
# for _, row in df_f.iterrows():
#     G.add_node(int(row["cluster index"]), **row.to_dict())
# for _, row in ef.iterrows():
#     G.add_edge(int(row["CLUSTERID1"]), int(row["CLUSTERID2"]),
#                cosine=row["Cosine"], delta_mz=row["DeltaMZ"])
#
# components = list(nx.connected_components(G))
# print(f"  Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")
# print(f"  Components: {len(components)}, Largest: {max(len(c) for c in components)}")
#
# # Only plot components with ≥2 nodes (singletons clutter the layout)
# multi_nodes = [n for c in components if len(c) >= 2 for n in c]
# H = G.subgraph(multi_nodes).copy()
# print(f"  Plotting {H.number_of_nodes()} nodes in multi-node components")
#
# pos = nx.spring_layout(H, seed=42, k=1.8, iterations=100)
#
# fig, ax = plt.subplots(figsize=(26, 20))
# ax.set_facecolor("#f5f5f5")
# fig.patch.set_facecolor("#f5f5f5")
#
# nx.draw_networkx_edges(H, pos, ax=ax, alpha=0.2, width=0.5, edge_color="#666666")
#
# NODE_R = 0.03   # node radius in data coordinates — increase if nodes overlap
#
# for node in H.nodes():
#     x, y = pos[node]
#     nd = G.nodes[node]
#     vals = [float(nd.get(c, 0) or 0) for c in ISOLATE_COLS]
#     total = sum(vals)
#
#     if total == 0:
#         circ = plt.Circle((x, y), NODE_R * 0.7, color="#cccccc", zorder=3,
#                            transform=ax.transData)
#         ax.add_patch(circ)
#     else:
#         fracs = [v / total for v in vals]
#         inset = ax.inset_axes(
#             [x - NODE_R, y - NODE_R, 2*NODE_R, 2*NODE_R],
#             transform=ax.transData
#         )
#         inset.pie(fracs, colors=COLORS, startangle=90,
#                   wedgeprops={"linewidth": 0.3, "edgecolor": "white"})
#         inset.set_aspect("equal")
#
#     # Label library-annotated nodes
#     name = nd.get("Compound_Name", "")
#     mqscore = nd.get("MQScore", 0)
#     if name and isinstance(name, str) and name not in ["", "nan"]:
#         short = name[:22] + "…" if len(name) > 22 else name
#         ax.text(x, y + NODE_R * 1.6, short, fontsize=5, ha="center",
#                 va="bottom", color="#1a1a1a", fontweight="bold", zorder=6)
#
# legend_patches = [
#     mpatches.Patch(color=COLORS[i], label=ISOLATE_LABELS[i])
#     for i in range(len(ISOLATE_COLS))
# ]
# ax.legend(handles=legend_patches, title="Isolate", loc="lower left",
#           fontsize=8, title_fontsize=9, framealpha=0.95, ncol=2,
#           frameon=True, edgecolor="#cccccc")
#
# ax.set_title(
#     "Classical Molecular Network — Marine Heterotrophic Bacterial Isolates\n"
#     "Broth-subtracted · Cosine ≥ 0.7 · Min 6 matched peaks",
#     fontsize=13, pad=16
# )
# ax.axis("off")
# plt.tight_layout()
# plt.savefig(OUT_DIR / "cmn_piechart_network.png", dpi=200, bbox_inches="tight")
# plt.savefig(OUT_DIR / "cmn_piechart_network.svg", bbox_inches="tight")
# plt.close()
# print("  Saved → analysis_output/cmn_piechart_network.png (.svg)")



# =============================================================================
# 10. EXPORT Cytoscape-ready files
# =============================================================================
# Node table for Cytoscape import
cyto_cols = ["cluster index", "number of spectra", "parent mass", "precursor charge",
             "RTMean", "origin", "match_quality", "flag_ion_mobility", "flag_high_ppm",
             "Compound_Name", "MQScore", "MZErrorPPM", "Adduct", "npclassifier_pathway",
             "npclassifier_class", "component"] + ISOLATE_COLS
cyto_cols = [c for c in cyto_cols if c in df_f.columns]
df_f[cyto_cols].to_csv(OUT_DIR / "cytoscape_node_table.tsv", sep="\t", index=False)
ef.to_csv(OUT_DIR / "cytoscape_edge_table.tsv", sep="\t", index=False)
print("Saved Cytoscape tables → analysis_output/cytoscape_node_table.tsv + edge_table.tsv")

print("\n" + "="*60)
print("ALL DONE — check analysis_output/")
print("="*60)