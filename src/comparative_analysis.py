# %% [markdown]
# The goal of this script is to perform quantiative analysis between arabidopsis thaliana varieties proteomes 

# %%
# =========================================================
# Isoform plots — setup & data loading
# =========================================================

from pathlib import Path
import pandas as pd

# ---------------------------------------------------------
# Base directory
# ---------------------------------------------------------

BASE_DIR = Path("/home/ggiordano/snap/main")

# ---------------------------------------------------------
# Input (normalized protein groups)
# ---------------------------------------------------------

DATA_DIR = BASE_DIR / "data" / "at_founders" / "comparative"
PG_NORM_DIR = BASE_DIR / "data" / "at_founders" / "orthology"/ "pg_normalized"

# ---------------------------------------------------------
# Output (figures)
# ---------------------------------------------------------

EXPORT_DIR = BASE_DIR / "export" / "at_founders" / "comparative"
EXPORT_DIR.mkdir(parents=True, exist_ok=True)

print("Loading normalized protein groups from:", PG_NORM_DIR)
print("Saving figures to:", EXPORT_DIR)


# %%
# =========================================================
# Isoform abundance analysis (publication-ready, clean + ticks)
# =========================================================

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as fm
import seaborn as sns
import pandas as pd
import numpy as np
import itertools
from pathlib import Path

# ---------------------------------------------------------
# Output directory
# ---------------------------------------------------------

EXPORT_DIR = Path("/home/ggiordano/snap/main/export/at_founders/comparative")
EXPORT_DIR.mkdir(parents=True, exist_ok=True)

print("Export directory:", EXPORT_DIR)

# ---------------------------------------------------------
# Fonts
# ---------------------------------------------------------

font_dir = Path.home() / ".local/share/fonts"

for font in font_dir.glob("Arial*.TTF"):
    fm.fontManager.addfont(str(font))

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 6,
    "axes.titlesize": 6,
    "axes.labelsize": 6,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "figure.dpi": 300
})

sns.set_style("white")

# ---------------------------------------------------------
# Figure size (IDENTICAL to your reference)
# ---------------------------------------------------------

mm_to_inch = 1 / 25.4
fig_width = 80 * mm_to_inch
fig_height = 50 * mm_to_inch

# =========================================================
# LOAD DATA
# =========================================================

pg_tables_norm = {}

pg_files = sorted(PG_NORM_DIR.glob("*_normalized.csv"))

if len(pg_files) == 0:
    raise FileNotFoundError(f"No normalized files found in {PG_NORM_DIR}")

for f in pg_files:
    acc = f.stem.split("_")[0]
    df = pd.read_csv(f)
    pg_tables_norm[acc] = df
    print(f"{acc}: {df.shape[0]} rows loaded")

print("Total accessions:", len(pg_tables_norm))

# =========================================================
# COMPUTE METRICS
# =========================================================

summary = []
all_data = []

for acc, df in pg_tables_norm.items():

    df = df.copy()

    proteins = df["Protein.Group"].astype(str)

    base = proteins.str.replace(r"_\d+$", "", regex=True)
    base = base.str.replace(r"\.\d+$", "", regex=True)

    df["base_gene"] = base

    ibaq_cols = [c for c in df.columns if "iBAQ" in c]

    if len(ibaq_cols) == 0:
        print(f"Warning: no iBAQ columns for {acc}")
        continue

    df["mean_abundance"] = np.log2(
        df[ibaq_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1) + 1
    )

    gene_groups = df.groupby("base_gene")

    sims = []
    ratios = []

    for gene, sub in gene_groups:

        if sub.shape[0] < 2:
            continue

        values = sub["mean_abundance"].values

        for i, j in itertools.combinations(range(len(values)), 2):

            x, y = values[i], values[j]

            diff = abs(x - y)
            sim = 1 / (1 + diff)

            sims.append(sim)
            ratios.append(diff)

            all_data.append({
                "accession": acc,
                "gene": gene,
                "similarity": sim,
                "log2_fc": diff
            })

    if len(sims) > 0:
        summary.append({
            "accession": acc,
            "mean_similarity": np.mean(sims),
            "median_similarity": np.median(sims),
            "mean_log2_fc": np.mean(ratios),
            "n_pairs": len(sims)
        })

    print(f"{acc}: {len(sims)} isoform pairs")

summary_df = pd.DataFrame(summary).sort_values("accession")
dist_df = pd.DataFrame(all_data)

# =========================================================
# COMMON AXIS SETTINGS (ticks visible but clean)
# =========================================================

def apply_ticks(ax):
    ax.tick_params(
        axis="both",
        which="major",
        direction="out",
        length=2,
        width=0.4
    )

# =========================================================
# PLOT 1 — BARPLOT
# =========================================================

fig, ax = plt.subplots(figsize=(fig_width, fig_height))

ax.bar(summary_df["accession"], summary_df["mean_similarity"],
       edgecolor="0.3", linewidth=0.15)

ax.set_ylabel("Mean isoform similarity")
ax.set_xlabel("Accession")
ax.set_ylim(0, 1)

ax.tick_params(axis="x", rotation=90)
apply_ticks(ax)

sns.despine()

plt.tight_layout()
plt.savefig(EXPORT_DIR / "plot1_mean_similarity.pdf")
plt.savefig(EXPORT_DIR / "plot1_mean_similarity.jpeg")
plt.close()

# =========================================================
# PLOT 2 — VIOLIN
# =========================================================

fig, ax = plt.subplots(figsize=(fig_width * 1.5, fig_height))

sns.violinplot(
    data=dist_df,
    x="accession",
    y="similarity",
    inner="quartile",
    cut=0,
    linewidth=0.3,
    ax=ax
)

ax.set_ylabel("Isoform similarity")
ax.set_xlabel("Accession")

ax.tick_params(axis="x", rotation=90)
apply_ticks(ax)

sns.despine()

plt.tight_layout()
plt.savefig(EXPORT_DIR / "plot2_violin_similarity.pdf")
plt.savefig(EXPORT_DIR / "plot2_violin_similarity.jpeg")
plt.close()

# =========================================================
# PLOT 3 — FRACTION DIVERGENT
# =========================================================

threshold = 0.5
dist_df["divergent"] = dist_df["similarity"] < threshold

div_df = dist_df.groupby("accession")["divergent"].mean().reset_index()

fig, ax = plt.subplots(figsize=(fig_width, fig_height))

ax.bar(div_df["accession"], div_df["divergent"],
       edgecolor="0.3", linewidth=0.15)

ax.set_ylabel("Fraction divergent isoforms")
ax.set_xlabel("Accession")
ax.set_ylim(0, 1)

ax.tick_params(axis="x", rotation=90)
apply_ticks(ax)

sns.despine()

plt.tight_layout()
plt.savefig(EXPORT_DIR / "plot3_divergent_fraction.pdf")
plt.savefig(EXPORT_DIR / "plot3_divergent_fraction.jpeg")
plt.close()

# =========================================================
# PLOT 4 — GLOBAL DISTRIBUTION
# =========================================================

fig, ax = plt.subplots(figsize=(fig_width, fig_height))

sns.histplot(dist_df["similarity"], bins=40, ax=ax)

ax.set_xlabel("Isoform similarity")
ax.set_ylabel("Count")

apply_ticks(ax)

sns.despine()

plt.tight_layout()
plt.savefig(EXPORT_DIR / "plot4_global_distribution.pdf")
plt.savefig(EXPORT_DIR / "plot4_global_distribution.jpeg")
plt.close()

# =========================================================
# PLOT 5 — LOG2 FC
# =========================================================

fig, ax = plt.subplots(figsize=(fig_width * 1.5, fig_height))

sns.boxplot(
    data=dist_df,
    x="accession",
    y="log2_fc",
    linewidth=0.3,
    fliersize=1,
    ax=ax
)

ax.set_ylabel("Log2 fold difference (isoforms)")
ax.set_xlabel("Accession")

ax.tick_params(axis="x", rotation=90)
apply_ticks(ax)

sns.despine()

plt.tight_layout()
plt.savefig(EXPORT_DIR / "plot5_log2fc.pdf")
plt.savefig(EXPORT_DIR / "plot5_log2fc.jpeg")
plt.close()

print("All plots saved in:", EXPORT_DIR)

# %%
for acc, df in pg_tables_norm.items():
    print("\n", acc)
    print(df.columns.tolist())
    break

# %%
# =========================================================
# Isoform abundance analysis (advanced, publication-ready)
# =========================================================

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as fm
import seaborn as sns
import pandas as pd
import numpy as np
import itertools
from pathlib import Path

# ---------------------------------------------------------
# Output directory
# ---------------------------------------------------------

EXPORT_DIR = Path("/home/ggiordano/snap/main/export/at_founders/comparative")
EXPORT_DIR.mkdir(parents=True, exist_ok=True)

print("Export directory:", EXPORT_DIR)

# ---------------------------------------------------------
# Fonts
# ---------------------------------------------------------

font_dir = Path.home() / ".local/share/fonts"

for font in font_dir.glob("Arial*.TTF"):
    fm.fontManager.addfont(str(font))

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 6,
    "axes.titlesize": 6,
    "axes.labelsize": 6,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "figure.dpi": 300
})

# ---------------------------------------------------------
# Figure size
# ---------------------------------------------------------

mm_to_inch = 1 / 25.4
fig_width = 80 * mm_to_inch
fig_height = 50 * mm_to_inch

# =========================================================
# COMPUTE METRICS
# =========================================================

summary = []
all_data = []

for acc, df in pg_tables_norm.items():

    df = df.copy()

    proteins = df["Protein.Group"].astype(str)
    base = proteins.str.replace(r"_\d+$", "", regex=True)
    base = base.str.replace(r"\.\d+$", "", regex=True)
    df["base_gene"] = base

    ibaq_cols = [c for c in df.columns if "iBAQ" in c]

    df["mean_abundance"] = np.log2(df[ibaq_cols].astype(float).mean(axis=1) + 1)

    gene_groups = df.groupby("base_gene")

    sims = []
    ratios = []

    for gene, sub in gene_groups:

        if sub.shape[0] < 2:
            continue

        values = sub["mean_abundance"].values

        for i, j in itertools.combinations(range(len(values)), 2):

            x, y = values[i], values[j]

            diff = abs(x - y)
            sim = 1 / (1 + diff)

            ratio = abs(x - y)

            sims.append(sim)
            ratios.append(ratio)

            all_data.append({
                "accession": acc,
                "gene": gene,
                "similarity": sim,
                "log2_fc": ratio
            })

    if len(sims) > 0:
        summary.append({
            "accession": acc,
            "mean_similarity": np.mean(sims),
            "median_similarity": np.median(sims),
            "mean_log2_fc": np.mean(ratios),
            "n_pairs": len(sims)
        })

    print(f"{acc}: {len(sims)} isoform pairs")

summary_df = pd.DataFrame(summary).sort_values("accession")
dist_df = pd.DataFrame(all_data)

# =========================================================
# COMMON TICK STYLE
# =========================================================

def apply_ticks(ax):
    # Tick styling
    ax.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=2,
        width=0.5,
        bottom=True,
        top=False,
        left=True,
        right=False
    )

    # Spine control: keep only bottom (x-axis frame)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_linewidth(0.5)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_linewidth(0.5)


# =========================================================
# PLOT 1 — BARPLOT
# =========================================================

fig, ax = plt.subplots(figsize=(fig_width, fig_height))

ax.bar(summary_df["accession"], summary_df["mean_similarity"],
       edgecolor="0.3", linewidth=0.15)

ax.set_ylabel("Mean isoform similarity")
ax.set_xlabel("Accession")
ax.set_ylim(0, 1)
ax.tick_params(axis="x", rotation=90)

apply_ticks(ax)

plt.tight_layout()
plt.savefig(EXPORT_DIR / "plot1_mean_similarity.pdf")
plt.savefig(EXPORT_DIR / "plot1_mean_similarity.jpeg", dpi=300)

plt.close()

# =========================================================
# PLOT 2 — VIOLIN
# =========================================================

fig, ax = plt.subplots(figsize=(fig_width*1.5, fig_height))

sns.violinplot(
    data=dist_df,
    x="accession",
    y="similarity",
    inner="quartile",
    cut=0,
    linewidth=0.3
)

ax.set_ylabel("Isoform similarity")
ax.set_xlabel("Accession")
ax.tick_params(axis="x", rotation=90)

apply_ticks(ax)

plt.tight_layout()
plt.savefig(EXPORT_DIR / "plot2_violin_similarity.pdf")
plt.savefig(EXPORT_DIR / "plot2_violin_similarity.jpeg", dpi=300)
plt.close()

# =========================================================
# PLOT 3 — FRACTION DIVERGENT
# =========================================================

threshold = 0.5
dist_df["divergent"] = dist_df["similarity"] < threshold

div_df = dist_df.groupby("accession")["divergent"].mean().reset_index()

fig, ax = plt.subplots(figsize=(fig_width, fig_height))

ax.bar(div_df["accession"], div_df["divergent"],
       edgecolor="0.3", linewidth=0.15)

ax.set_ylabel("Fraction divergent isoforms")
ax.set_xlabel("Accession")
ax.set_ylim(0, 1)
ax.tick_params(axis="x", rotation=90)

apply_ticks(ax)

plt.tight_layout()
plt.savefig(EXPORT_DIR / "plot3_divergent_fraction.pdf")
plt.savefig(EXPORT_DIR / "plot3_divergent_fraction.jpeg", dpi=300)
plt.close()

# =========================================================
# PLOT 4 — GLOBAL DISTRIBUTION
# =========================================================

fig, ax = plt.subplots(figsize=(fig_width, fig_height))

sns.histplot(dist_df["similarity"], bins=40)

ax.set_xlabel("Isoform similarity")
ax.set_ylabel("Count")

apply_ticks(ax)

plt.tight_layout()
plt.savefig(EXPORT_DIR / "plot4_global_distribution.pdf")
plt.savefig(EXPORT_DIR / "plot4_global_distribution.jpeg", dpi=300)
plt.close()

# =========================================================
# PLOT 5 — LOG2 FOLD DIFFERENCE
# =========================================================

fig, ax = plt.subplots(figsize=(fig_width*1.5, fig_height))

sns.boxplot(
    data=dist_df,
    x="accession",
    y="log2_fc",
    linewidth=0.3,
    fliersize=1
)

ax.set_ylabel("Log2 fold difference (isoforms)")
ax.set_xlabel("Accession")
ax.tick_params(axis="x", rotation=90)

apply_ticks(ax)

plt.tight_layout()
plt.savefig(EXPORT_DIR / "plot5_log2fc.pdf")
plt.savefig(EXPORT_DIR / "plot5_log2fc.jpeg", dpi=300)

plt.close()

print("All plots saved in:", EXPORT_DIR)

# %%
# =========================================================
# Core / Shell / Cloud per accession (clean + fixed)
# =========================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as fm
from pathlib import Path

# ---------------------------------------------------------
# INPUT / OUTPUT
# ---------------------------------------------------------

INPUT_FILE = "/home/ggiordano/snap/main/data/at_founders/orthology/orthology_with_identified_proteins/orthology_matrix_final.tsv"
EXPORT_DIR = Path("/home/ggiordano/snap/main/export/at_founders/comparative")
EXPORT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------
# FONTS
# ---------------------------------------------------------

font_dir = Path.home() / ".local/share/fonts"
for font in font_dir.glob("Arial*.TTF"):
    fm.fontManager.addfont(str(font))

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

# ---------------------------------------------------------
# FIGURE SIZE
# ---------------------------------------------------------

mm_to_inch = 1 / 25.4
fig_width = 180 * mm_to_inch
fig_height = 80 * mm_to_inch

# =========================================================
# LOAD DATA
# =========================================================

df = pd.read_csv(INPUT_FILE, sep="\t")

# =========================================================
# DETECT ACCESSIONS
# =========================================================

accessions = sorted(set(
    col.split("_")[0]
    for col in df.columns
    if "_iBAQ_" in col
))

print("Detected accessions:", accessions)

# =========================================================
# PARSE iBAQ
# =========================================================

def parse_ibaq(cell):
    if pd.isna(cell):
        return 0.0
    if isinstance(cell, (int, float)):
        return float(cell)
    vals = [v.strip() for v in str(cell).split(",")]
    nums = []
    for v in vals:
        try:
            nums.append(float(v))
        except:
            continue
    return max(nums) if nums else 0.0

# =========================================================
# BUILD PRESENCE MATRIX
# =========================================================

presence = {}

for acc in accessions:
    ibaq_cols = [c for c in df.columns if c.startswith(f"{acc}_iBAQ_")]
    parsed = df[ibaq_cols].apply(lambda col: col.map(parse_ibaq))
    presence[acc] = parsed.max(axis=1) > 0

presence_df = pd.DataFrame(presence)
presence_df.index = df["Orthogroup"]

# =========================================================
# CLASSIFICATION
# =========================================================

counts = presence_df.sum(axis=1)
n_acc = len(accessions)

core_mask  = counts == n_acc
shell_mask = (counts > 1) & (counts < n_acc)
cloud_mask = counts == 1

# =========================================================
# SUMMARY TABLE
# =========================================================

summary = []

for acc in accessions:
    summary.append({
        "accession": acc,
        "Core":  (presence_df[acc] & core_mask).sum(),
        "Shell": (presence_df[acc] & shell_mask).sum(),
        "Cloud": (presence_df[acc] & cloud_mask).sum()
    })

summary_df = pd.DataFrame(summary)

# =========================================================
# PLOT — ABSOLUTE COUNTS WITH PROPER BROKEN AXIS
# =========================================================

fig, (ax_top, ax_bottom) = plt.subplots(
    2, 1, sharex=True,
    gridspec_kw={'height_ratios': [1.2, 1]},
    figsize=(fig_width, fig_height)
)

x = np.arange(len(summary_df))

core_vals = summary_df["Core"].values
shell_vals = summary_df["Shell"].values
cloud_vals = summary_df["Cloud"].values
total_vals = core_vals + shell_vals + cloud_vals

colors = ["#2c3e50", "#3498db", "#e74c3c"]

# ---------------------------------------------------------
# BARS (NO WHITE GAPS)
# ---------------------------------------------------------

bars_core = ax_bottom.bar(
    x, core_vals,
    color=colors[0],
    edgecolor="none",
    linewidth=0
)

bars_shell = ax_bottom.bar(
    x, shell_vals,
    bottom=core_vals,
    color=colors[1],
    edgecolor="none",
    linewidth=0
)

bars_cloud = ax_bottom.bar(
    x, cloud_vals,
    bottom=core_vals + shell_vals,
    color=colors[2],
    edgecolor="none",
    linewidth=0
)

# replicate on top axis (no handles needed)
for ax in [ax_top]:
    ax.bar(x, core_vals, color=colors[0], edgecolor="none", linewidth=0)
    ax.bar(x, shell_vals, bottom=core_vals, color=colors[1], edgecolor="none", linewidth=0)
    ax.bar(x, cloud_vals, bottom=core_vals + shell_vals, color=colors[2], edgecolor="none", linewidth=0)

# ---------------------------------------------------------
# STRONG CUT (FOCUS ON CLOUD)
# ---------------------------------------------------------

cut = 7000        # where compression happens
gap = 200         # small gap only (IMPORTANT)

ax_bottom.set_ylim(0, cut - gap)
ax_top.set_ylim(cut + gap, total_vals.max() * 1.02)

# ---------------------------------------------------------
# STYLE — CLEAN (NO DIVIDING LINES)
# ---------------------------------------------------------

# remove ALL unnecessary spines
for ax in [ax_top, ax_bottom]:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

# keep only left + bottom where needed
ax_bottom.spines['left'].set_linewidth(0.4)
ax_bottom.spines['bottom'].set_linewidth(0.4)

ax_top.spines['left'].set_linewidth(0.4)
ax_top.spines['bottom'].set_visible(False)   # ← IMPORTANT (removes dividing line)

# ---------------------------------------------------------
# TICKS
# ---------------------------------------------------------

for ax in [ax_top, ax_bottom]:
    ax.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=2,
        width=0.4,
        bottom=True,
        top=False,
        left=True,
        right=False
    )

# remove top panel x ticks completely (cleaner)
ax_top.tick_params(bottom=False)

# ticks
for ax in [ax_top, ax_bottom]:
    ax.tick_params(width=0.4, length=2)

# ---------------------------------------------------------
# BREAK MARKS
# ---------------------------------------------------------

d = .005
kwargs = dict(transform=ax_top.transAxes, color='k', clip_on=False, lw=0.5)

# LEFT ONLY
ax_top.plot((-d, +d), (-d, +d), **kwargs)

kwargs.update(transform=ax_bottom.transAxes)
ax_bottom.plot((-d, +d), (1 - d, 1 + d), **kwargs)

# ---------------------------------------------------------
# LABELS
# ---------------------------------------------------------

ax_bottom.set_xticks(x)
ax_bottom.set_xticklabels(summary_df["accession"])

fig.text(0.04, 0.5, 'Number of orthogroups', va='center', rotation='vertical')

# ---------------------------------------------------------
# LEGEND (CORRECT)
# ---------------------------------------------------------

fig.legend(
    [bars_core, bars_shell, bars_cloud],
    ["Core", "Shell", "Cloud"],
    loc="upper right",
    bbox_to_anchor=(0.98, 0.98),
    frameon=False
)

plt.subplots_adjust(hspace=0.05)

plt.savefig(EXPORT_DIR / "core_shell_cloud_broken_axis.pdf", bbox_inches="tight")
plt.savefig(EXPORT_DIR / "core_shell_cloud_broken_axis.jpeg", bbox_inches="tight")

print("Plot saved")

plt.show()

# =========================================================
# COMPACT PROTEIN EXPORT
# =========================================================

COMPACT_EXPORT = Path("/home/ggiordano/snap/main/data/at_founders/comparative")
COMPACT_EXPORT.mkdir(parents=True, exist_ok=True)

protein_cols = {
    acc: [c for c in df.columns if c.startswith(f"{acc}_") and "_iBAQ_" not in c]
    for acc in accessions
}

rows = []

for acc in accessions:

    acc_present = presence_df[acc]

    class_masks = {
        "core":  acc_present & core_mask,
        "shell": acc_present & shell_mask,
        "cloud": acc_present & cloud_mask
    }

    cols = protein_cols[acc]

    for cls, mask in class_masks.items():

        subset = df.loc[mask.values, cols]

        proteins = set()

        for col in subset.columns:
            for val in subset[col].dropna():
                parts = [p.strip() for p in str(val).split(",")]
                for p in parts:
                    if p and p != "nan":
                        proteins.add(p)

        rows.append({
            "accession": acc,
            "class": cls,
            "n_proteins": len(proteins),
            "proteins": ";".join(sorted(proteins))
        })

protein_table = pd.DataFrame(rows)

out_file = COMPACT_EXPORT / "core_shell_cloud_proteins_compact.tsv"
protein_table.to_csv(out_file, sep="\t", index=False)

print("Saved:", out_file)

# FIXED PRINT
print(summary_df[["accession", "Cloud"]])

# %%
# =========================================================
# DIFFERENTIAL ABUNDANCE (ORTHOGROUP LEVEL — LIMMA-LIKE)
# =========================================================

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import seaborn as sns

print("\n--- Differential abundance analysis (limma-like) ---")

# ---------------------------------------------------------
# BUILD ORTHOGROUP MATRIX
# ---------------------------------------------------------

long_data = []

for acc, df in pg_tables_norm.items():

    df = df.copy()

    proteins = df["Protein.Group"].astype(str)
    base = proteins.str.replace(r"_\d+$", "", regex=True)
    base = base.str.replace(r"\.\d+$", "", regex=True)
    df["orthogroup"] = base

    ibaq_cols = [c for c in df.columns if "iBAQ" in c]

    for col in ibaq_cols:

        tmp = pd.DataFrame({
            "orthogroup": df["orthogroup"],
            "value": pd.to_numeric(df[col], errors="coerce"),
            "sample": f"{acc}_{col}",
            "accession": acc
        })

        long_data.append(tmp)

long_df = pd.concat(long_data)

# sum isoforms → orthogroup
long_df = long_df.groupby(["orthogroup", "sample", "accession"])["value"].sum().reset_index()

matrix = long_df.pivot_table(
    index="orthogroup",
    columns="sample",
    values="value"
)

metadata = long_df[["sample", "accession"]].drop_duplicates().set_index("sample")

# ---------------------------------------------------------
# LOG2 + IMPUTATION
# ---------------------------------------------------------

matrix = np.log2(matrix + 1)

global_min = np.nanpercentile(matrix.values, 1)
matrix = matrix.fillna(global_min)

# ---------------------------------------------------------
# DIFFERENTIAL TEST (ONE-VS-ALL)
# ---------------------------------------------------------

results = []

for og in matrix.index:

    values = matrix.loc[og]
    values = values.dropna()

    for acc in metadata["accession"].unique():

        group_mask = metadata.loc[values.index, "accession"] == acc

        group_vals = values[group_mask]
        other_vals = values[~group_mask]

        # skip if too few points
        if len(group_vals) < 2 or len(other_vals) < 2:
            continue

        # log2FC
        log2fc = group_vals.mean() - other_vals.mean()

        # t-test
        stat, pval = ttest_ind(group_vals, other_vals, equal_var=False)

        results.append({
            "orthogroup": og,
            "accession": acc,
            "log2FC": log2fc,
            "pval": pval
        })

res_df = pd.DataFrame(results)

# ---------------------------------------------------------
# MULTIPLE TESTING
# ---------------------------------------------------------

res_df["adj_pval"] = multipletests(res_df["pval"], method="fdr_bh")[1]

# ---------------------------------------------------------
# SIGNIFICANT
# ---------------------------------------------------------

logfc_threshold = 1
alpha = 0.05

sig_df = res_df[
    (res_df["adj_pval"] < alpha) &
    (np.abs(res_df["log2FC"]) > logfc_threshold)
]

# ---------------------------------------------------------
# SAVE
# ---------------------------------------------------------

OUT_DIR = EXPORT_DIR / "differential_analysis"
OUT_DIR.mkdir(exist_ok=True)

res_df.to_csv(OUT_DIR / "all_results.tsv", sep="\t", index=False)
sig_df.to_csv(OUT_DIR / "significant_orthogroups.tsv", sep="\t", index=False)

print("Saved results to:", OUT_DIR)

# ---------------------------------------------------------
# SUMMARY
# ---------------------------------------------------------

summary_hits = sig_df.groupby("accession").size().reset_index(name="n_hits")
print("\nSignificant orthogroups per accession:")
print(summary_hits)

# =========================================================
# VOLCANO PLOTS
# =========================================================

print("\n--- Generating volcano plots ---")

for acc in res_df["accession"].unique():

    sub = res_df[res_df["accession"] == acc].copy()

    sub["neglog10"] = -np.log10(sub["adj_pval"])

    # color logic
    sub["color"] = "lightblue"
    sub.loc[
        (sub["adj_pval"] < alpha) & (np.abs(sub["log2FC"]) > logfc_threshold),
        "color"
    ] = "red"

    plt.figure(figsize=(4, 3))

    plt.scatter(
        sub["log2FC"],
        sub["neglog10"],
        c=sub["color"],
        s=5,
        alpha=0.7
    )

    plt.axvline(logfc_threshold, linestyle="--", linewidth=0.5)
    plt.axvline(-logfc_threshold, linestyle="--", linewidth=0.5)
    plt.axhline(-np.log10(alpha), linestyle="--", linewidth=0.5)

    plt.xlabel("log2 fold change")
    plt.ylabel("-log10 adj p-value")
    plt.title(acc)

    plt.tight_layout()

    plt.savefig(OUT_DIR / f"volcano_{acc}.pdf")
    plt.savefig(OUT_DIR / f"volcano_{acc}.jpeg", dpi=300)

    plt.close()

print("Volcano plots saved.")

# %%
print("Design shape:", design.shape)
print("Matrix shape:", matrix_imputed.shape)
print(res_df.head())
print(res_df.columns)


