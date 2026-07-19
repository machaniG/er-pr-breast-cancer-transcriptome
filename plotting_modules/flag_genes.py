from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def load_deseq_table(csv_path, gene_col="gene_id"):
    df = pd.read_csv(csv_path, index_col=0).copy()
    if df.index.name is None:
        df.index.name = gene_col
    df = df.rename_axis(gene_col).reset_index()
    df["padj"] = df["padj"].fillna(1.0)
    df["log2FoldChange"] = df["log2FoldChange"].fillna(0.0)
    return df


def merge_wt_mutant_tables(
    wt_csv,
    mutant_csv,
    gene_col="gene_id",
    wt_prefix="wt",
    mutant_prefix="mutant",
):
    wt = load_deseq_table(wt_csv, gene_col=gene_col).rename(
        columns={
            "padj": f"{wt_prefix}_padj",
            "log2FoldChange": f"{wt_prefix}_log2FoldChange",
        }
    )
    mutant = load_deseq_table(mutant_csv, gene_col=gene_col).rename(
        columns={
            "padj": f"{mutant_prefix}_padj",
            "log2FoldChange": f"{mutant_prefix}_log2FoldChange",
        }
    )

    merged = mutant.merge(
        wt[[gene_col, f"{wt_prefix}_padj", f"{wt_prefix}_log2FoldChange"]],
        on=gene_col,
        how="left",
    )

    return merged


def add_mutant_specific_flag(
    df,
    wt_padj_col="wt_padj",
    wt_log2fc_col="wt_log2FoldChange",
    mutant_padj_col="mutant_padj",
    mutant_log2fc_col="mutant_log2FoldChange",
    padj_threshold=0.05,
    log2fc_threshold=1.0,
    wt_effect_threshold=0.5,
):
    df = df.copy()
    df["mutant_sig_down"] = (
        (df[mutant_padj_col] < padj_threshold)
        & (df[mutant_log2fc_col] <= -log2fc_threshold)
    )
    df["wt_untouched"] = (
        (df[wt_padj_col].isna())
        | (
            (df[wt_padj_col] >= padj_threshold)
            & (df[wt_log2fc_col].abs() < wt_effect_threshold)
        )
    )
    df["mutant_specific_target"] = df["mutant_sig_down"] & df["wt_untouched"]
    df["minus_log10_padj"] = -np.log10(df[mutant_padj_col].clip(lower=1e-300))
    return df


def plot_mutant_volcano_with_highlights(
    merged_df,
    out_png,
    gene_col="gene_id",
    mutant_log2fc_col="mutant_log2FoldChange",
    mutant_padj_col="mutant_padj",
    padj_threshold=0.05,
    log2fc_threshold=1.0,
    highlight_col="mutant_specific_target",
    label_top_n=10,
    title="shPR knockdown effect in Y537S mutant",
    figsize=(8, 7),
    point_size=14,
):
    df = merged_df.copy()

    sig = (df[mutant_padj_col] < padj_threshold) & (
        df[mutant_log2fc_col].abs() >= log2fc_threshold
    )
    up = sig & (df[mutant_log2fc_col] > 0)
    down = sig & (df[mutant_log2fc_col] < 0)
    ns = ~sig
    hi = df[highlight_col].fillna(False)

    fig, ax = plt.subplots(figsize=figsize)

    ax.scatter(
        df.loc[ns & ~hi, mutant_log2fc_col],
        df.loc[ns & ~hi, "minus_log10_padj"],
        s=point_size,
        c="#C0C0C0",
        alpha=0.5,
        linewidths=0,
        label="Not significant",
    )
    ax.scatter(
        df.loc[down & ~hi, mutant_log2fc_col],
        df.loc[down & ~hi, "minus_log10_padj"],
        s=point_size,
        c="#4C78A8",
        alpha=0.8,
        linewidths=0,
        label="Down",
    )
    ax.scatter(
        df.loc[up & ~hi, mutant_log2fc_col],
        df.loc[up & ~hi, "minus_log10_padj"],
        s=point_size,
        c="#D62728",
        alpha=0.8,
        linewidths=0,
        label="Up",
    )
    ax.scatter(
        df.loc[hi, mutant_log2fc_col],
        df.loc[hi, "minus_log10_padj"],
        s=point_size * 1.8,
        c="#8A2BE2",
        alpha=0.95,
        linewidths=0.5,
        edgecolors="black",
        label="Mutant-specific candidate",
    )

    ax.axvline(-log2fc_threshold, color="black", linestyle="--", linewidth=1)
    ax.axvline(log2fc_threshold, color="black", linestyle="--", linewidth=1)
    ax.axhline(-np.log10(padj_threshold), color="black", linestyle="--", linewidth=1)

    top = df.loc[hi].sort_values(mutant_padj_col).head(label_top_n)
    for _, row in top.iterrows():
        ax.text(
            row[mutant_log2fc_col],
            row["minus_log10_padj"],
            str(row[gene_col]),
            fontsize=8,
            ha="left",
            va="bottom",
        )

    ax.set_title(title)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(fig)

    return out_png


def save_candidate_table(merged_df, out_csv, gene_col="gene_id"):
    cols = [
        gene_col,
        "mutant_log2FoldChange",
        "mutant_padj",
        "wt_log2FoldChange",
        "wt_padj",
        "mutant_specific_target",
    ]
    out = merged_df.loc[merged_df["mutant_specific_target"], cols].copy()
    out = out.sort_values("mutant_padj")
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_csv, index=False)
    return out_csv