# =============================================================================
# 02_broth_shared_analysis.py
# MicroChemEco Lab | R. Wang
#
# Complementary to 00_gnps_cmn_analysis.py — instead of removing broth-related
# nodes, this script looks ONLY at them: broth-only nodes + nodes shared
# between broth and isolates. Produces one summary figure of library matches
# (and m/z + RT for unidentified hits) plus the underlying detail tables.
#
# Input: nf_output/networking/clustersummary_with_network.tsv
# =============================================================================

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

REPO_ROOT = Path(__file__).resolve().parents[1]
NETWORK_DIR = REPO_ROOT / "data" / "MCE-GNPS" / "nf_output" / "networking"
OUT_DIR = REPO_ROOT / "results" / "msms_analysis_output"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MAIN_FILE = NETWORK_DIR / "clustersummary_with_network.tsv"

BROTH_COL = "ATTRIBUTE_PERFILEGROUPING:GNPSGROUP:bkgrnd-marine-broth.mzML"
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

# ── 1. LOAD ───────────────────────────────────────────────────────────────────
print("Loading main node table...")
df = pd.read_csv(MAIN_FILE, sep="\t", low_memory=False)
print(f"  {len(df)} nodes total")

# ── 2. SPLIT INTO BROTH-ONLY / SHARED ────────────────────────────────────────
broth_present = df[BROTH_COL] > 0
isolate_present = df[ISOLATE_COLS].sum(axis=1) > 0

broth_only = df[broth_present & ~isolate_present].copy()
shared = df[broth_present & isolate_present].copy()

print(f"  Broth-only nodes:            {len(broth_only)}")
print(f"  Shared (broth + isolate):    {len(shared)}")

# ── 3. CLASSIFY LIBRARY MATCHES (same flagging rule as 01_eda_data_cleaning.py) ──
PPM_FLAG_LIMIT = 300

def classify(d):
    has_match = d["Compound_Name"].notna() & (d["Compound_Name"] != "")
    flag_ion_mobility = has_match & d.get("LibraryName", pd.Series("", index=d.index)).str.contains("ION-MOBILITY", na=False)
    flag_high_ppm = has_match & (d["MZErrorPPM"].abs() > PPM_FLAG_LIMIT)
    flagged = flag_ion_mobility | flag_high_ppm
    d = d.copy()
    d["flag_ion_mobility"] = flag_ion_mobility
    d["flag_high_ppm"] = flag_high_ppm
    d["match_quality"] = "UNANNOTATED"
    d.loc[has_match & ~flagged, "match_quality"] = "VALID"
    d.loc[flagged, "match_quality"] = "FLAGGED"
    return d

broth_only = classify(broth_only)
shared = classify(shared)

print("\nBroth-only match status:")
print(broth_only["match_quality"].value_counts().to_string())
print(f"  (of which ion-mobility: {broth_only['flag_ion_mobility'].sum()}, high ppm: {broth_only['flag_high_ppm'].sum()})")
print("\nShared match status:")
print(shared["match_quality"].value_counts().to_string())
print(f"  (of which ion-mobility: {shared['flag_ion_mobility'].sum()}, high ppm: {shared['flag_high_ppm'].sum()})")

# ── 4. SAVE DETAIL TABLES ─────────────────────────────────────────────────────
# RTMean from GNPS2 is reported in SECONDS (confirmed against RTINSECONDS in
# clustering/specs_ms.mgf), despite the dashboard column header saying "(Min)".
# RTMean_min below is the corrected value for downstream use.
broth_only["RTMean_min"] = broth_only["RTMean"] / 60
shared["RTMean_min"] = shared["RTMean"] / 60

detail_cols = [c for c in [
    "cluster index", "number of spectra", "parent mass", "RTMean_min",
    "Compound_Name", "MQScore", "MZErrorPPM", "SharedPeaks", "Adduct",
    "LibraryName", "npclassifier_pathway", "npclassifier_class",
    "match_quality", "flag_ion_mobility", "flag_high_ppm",
] if c in df.columns or c == "RTMean_min" or c in ("match_quality", "flag_ion_mobility", "flag_high_ppm")]

broth_only[detail_cols].to_csv(OUT_DIR / "broth_only_nodes.tsv", sep="\t", index=False)
shared[detail_cols].to_csv(OUT_DIR / "shared_broth_isolate_nodes.tsv", sep="\t", index=False)
print("\nSaved -> analysis_output/broth_only_nodes.tsv")
print("Saved -> analysis_output/shared_broth_isolate_nodes.tsv")

# ── 5. SUMMARY FIGURE (shared broth+isolate nodes) ───────────────────────────
valid_shared = shared[shared["match_quality"] == "VALID"]
n_valid = len(valid_shared)
n_flagged = (shared["match_quality"] == "FLAGGED").sum()
n_unidentified = (shared["match_quality"] == "UNANNOTATED").sum()

top_matches = (
    valid_shared["Compound_Name"]
    .value_counts()
    .head(10)
    .sort_values()
)

fig, axes = plt.subplots(1, 2, figsize=(15, 6))

# Panel A — top library matches among shared nodes
ax = axes[0]
if len(top_matches):
    top_matches.plot(kind="barh", ax=ax, color="#1f78b4")
    ax.set_xlabel("Number of nodes")
    ax.set_title(
        f"Top library matches — shared broth+isolate nodes\n"
        f"(VALID={n_valid}, FLAGGED={n_flagged}, unidentified={n_unidentified})"
    )
else:
    ax.text(0.5, 0.5, "No valid library matches", ha="center", va="center")
    ax.set_axis_off()

# Panel B — m/z vs RT for all shared nodes, colored by match status
# NOTE: RTMean from GNPS2 is in SECONDS despite the dashboard column header
# saying "(Min)" -- confirmed against RTINSECONDS in clustering/specs_ms.mgf.
ax = axes[1]
status_colors = {"VALID": "#2196F3", "FLAGGED": "#FF9800", "UNANNOTATED": "#BDBDBD"}
for status, color in status_colors.items():
    sub = shared[shared["match_quality"] == status]
    ax.scatter(sub["RTMean"] / 60, sub["parent mass"], s=22, alpha=0.7,
               color=color, label=f"{status} (n={len(sub)})", edgecolors="white", linewidths=0.3)
ax.set_xlabel("Retention time, RTMean (min)")
ax.set_ylabel("Parent mass (m/z)")
ax.set_title("m/z vs RT — shared broth+isolate nodes")
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig(OUT_DIR / "broth_shared_library_hits.png", dpi=150, bbox_inches="tight")
plt.close()
print("\nSaved -> analysis_output/broth_shared_library_hits.png")

print("\n" + "=" * 60)
print("ALL DONE — check analysis_output/")
print("=" * 60)
