#!/usr/bin/env python3
"""Generate a clean, 3-color volcano plot from DESeq2-style results.

Usage as a script:
    python plot_volcano.py --input deseq_results/treatment_shPR_vs_shGFP.up.csv \
        --output volcano.png

Or import `plot_volcano(res_df)` directly in a notebook.
"""

from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_volcano(
    res_df: pd.DataFrame,
    log2fc_col: str = "log2FoldChange",
    padj_col: str = "padj",
    label_col: str | None = "symbol",
    log2fc_threshold: float = 1.0,
    padj_threshold: float = 0.05,
    top_n_labels: int = 15,
    clip_log2fc: tuple[float, float] | None = None,
    title: str = "Volcano Plot of Differential Expression",
    figsize: tuple[float, float] = (8, 6.5),
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Plot a 3-color (up/down/not-significant) volcano plot.

    Points with padj == 0 are clipped to the smallest nonzero padj in the
    dataset (divided by 10) instead of producing -log10(0) = inf, which is
    what was pushing points to the extreme top of the original plot.

    `clip_log2fc`, e.g. (-25, 25), caps extreme log2FC outliers to that range
    for plotting so a single far-out point (like one gene at -40 among others
    capped around -25) doesn't compress the rest of the plot. Clipped points
    are drawn as hollow triangles at the boundary, distinct from real values,
    so they're never mistaken for a gene that actually sits at that x-position.
    """
    df = res_df.dropna(subset=[log2fc_col, padj_col]).copy()

    # Avoid -log10(0) = inf: clip zero p-values to just below the smallest observed nonzero value
    nonzero_min = df.loc[df[padj_col] > 0, padj_col].min()
    floor = nonzero_min / 10 if pd.notna(nonzero_min) else 1e-300
    df["padj_plot"] = df[padj_col].clip(lower=floor)
    df["neg_log10_padj"] = -np.log10(df["padj_plot"])

    # Avoid one extreme log2FC outlier compressing the x-axis for everything else
    df["log2fc_plot"] = df[log2fc_col]
    df["is_clipped"] = False
    if clip_log2fc is not None:
        lo, hi = clip_log2fc
        df["is_clipped"] = (df[log2fc_col] < lo) | (df[log2fc_col] > hi)
        df["log2fc_plot"] = df[log2fc_col].clip(lower=lo, upper=hi)

    is_sig = df[padj_col] < padj_threshold
    is_up = is_sig & (df[log2fc_col] > log2fc_threshold)
    is_down = is_sig & (df[log2fc_col] < -log2fc_threshold)
    is_ns = ~(is_up | is_down)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize, dpi=150)

    plot_col = "log2fc_plot"
    not_clipped = ~df["is_clipped"]

    # Order matters: plot NS first so colored points sit on top
    ax.scatter(
        df.loc[is_ns & not_clipped, plot_col], df.loc[is_ns & not_clipped, "neg_log10_padj"],
        c="#B0B0B0", s=12, alpha=0.4, linewidths=0, label="Not significant", rasterized=True,
    )
    ax.scatter(
        df.loc[is_down & not_clipped, plot_col], df.loc[is_down & not_clipped, "neg_log10_padj"],
        c="#2166AC", s=16, alpha=0.75, linewidths=0, label=f"Down ({is_down.sum()})", rasterized=True,
    )
    ax.scatter(
        df.loc[is_up & not_clipped, plot_col], df.loc[is_up & not_clipped, "neg_log10_padj"],
        c="#B2182B", s=16, alpha=0.75, linewidths=0, label=f"Up ({is_up.sum()})", rasterized=True,
    )

    # Clipped outliers: hollow triangles pointing outward, at the axis boundary
    if df["is_clipped"].any():
        clipped = df[df["is_clipped"]]
        for _, row in clipped.iterrows():
            color = "#B2182B" if row[log2fc_col] > 0 else "#2166AC"
            marker = ">" if row[log2fc_col] > 0 else "<"
            ax.scatter(
                row[plot_col], row["neg_log10_padj"],
                marker=marker, s=70, facecolors="white", edgecolors=color, linewidths=1.5, zorder=5,
            )
        ax.scatter([], [], marker=">", facecolors="white", edgecolors="black", linewidths=1.5,
                   label=f"Clipped outlier(s) ({df['is_clipped'].sum()})")

    ax.axhline(-np.log10(padj_threshold), color="#444444", linestyle="--", linewidth=0.8, zorder=0)
    ax.axvline(log2fc_threshold, color="#444444", linestyle="--", linewidth=0.8, zorder=0)
    ax.axvline(-log2fc_threshold, color="#444444", linestyle="--", linewidth=0.8, zorder=0)

    # Label the most extreme significant genes (by combined effect size)
    resolved_label_col = label_col
    if label_col is not None and label_col not in df.columns:
        fallback_names = ["symbol", "gene_symbol", "Gene", "gene_name", "GeneSymbol",
                           "external_gene_name", "hgnc_symbol", "Symbol", "gene"]
        resolved_label_col = next((c for c in fallback_names if c in df.columns), None)
        if resolved_label_col is not None:
            print(f"[plot_volcano] '{label_col}' column not found; using '{resolved_label_col}' for gene labels instead.")
        else:
            print(
                f"[plot_volcano] WARNING: no label column found (tried '{label_col}' and "
                f"{fallback_names}). No gene labels will be drawn. Available columns: {list(df.columns)}. "
                f"Pass label_col='<your column>' explicitly to fix this."
            )

    if resolved_label_col is not None and top_n_labels > 0:
        sig_df = df.loc[is_up | is_down].copy()
        if sig_df[resolved_label_col].isna().all():
            print(f"[plot_volcano] WARNING: label column '{resolved_label_col}' is entirely empty for significant genes; no labels drawn.")
        sig_df["score"] = sig_df[log2fc_col].abs() * sig_df["neg_log10_padj"]
        top_genes = sig_df.nlargest(top_n_labels, "score")
        for _, row in top_genes.iterrows():
            gene_label = row[resolved_label_col]
            if pd.isna(gene_label):
                continue
            ax.annotate(
                str(gene_label),
                (row[plot_col], row["neg_log10_padj"]),
                fontsize=7, color="black",
                xytext=(3, 3), textcoords="offset points",
            )

    ax.set_title(title, fontsize=13, weight="bold")
    ax.set_xlabel("log$_2$ Fold Change", fontsize=11)
    ax.set_ylabel("-log$_{10}$ Adjusted P-value", fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, loc="upper right", fontsize=9, markerscale=1.5)
    ax.margins(x=0.05, y=0.05)
    plt.tight_layout()
    return ax


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate a 3-color volcano plot from a DESeq2 results CSV")
    parser.add_argument("--input", required=True, help="CSV with log2FoldChange and padj columns")
    parser.add_argument("--output", default="volcano.png", help="Output image path")
    parser.add_argument("--log2fc-col", default="log2FoldChange")
    parser.add_argument("--padj-col", default="padj")
    parser.add_argument("--label-col", default="symbol", help="Column for gene labels (set to '' to disable)")
    parser.add_argument("--log2fc-threshold", type=float, default=1.0)
    parser.add_argument("--padj-threshold", type=float, default=0.05)
    parser.add_argument("--top-n-labels", type=int, default=15)
    parser.add_argument("--clip-log2fc-min", type=float, default=None, help="Clip extreme negative log2FC outliers to this value")
    parser.add_argument("--clip-log2fc-max", type=float, default=None, help="Clip extreme positive log2FC outliers to this value")
    parser.add_argument("--title", default="Volcano Plot of Differential Expression")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    res_df = pd.read_csv(args.input)
    clip_log2fc = None
    if args.clip_log2fc_min is not None or args.clip_log2fc_max is not None:
        lo = args.clip_log2fc_min if args.clip_log2fc_min is not None else res_df[args.log2fc_col].min()
        hi = args.clip_log2fc_max if args.clip_log2fc_max is not None else res_df[args.log2fc_col].max()
        clip_log2fc = (lo, hi)

    plot_volcano(
        res_df,
        log2fc_col=args.log2fc_col,
        padj_col=args.padj_col,
        label_col=args.label_col or None,
        log2fc_threshold=args.log2fc_threshold,
        padj_threshold=args.padj_threshold,
        top_n_labels=args.top_n_labels,
        clip_log2fc=clip_log2fc,
        title=args.title,
    )
    plt.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"Saved volcano plot to {args.output}")


if __name__ == "__main__":
    main()
