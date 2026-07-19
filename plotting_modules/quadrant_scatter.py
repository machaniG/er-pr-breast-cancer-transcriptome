import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_quadrant_scatter(
    df_x,
    df_y,
    label_x="contrast_x",
    label_y="contrast_y",
    gene_col=None,
    lfc_col="log2FoldChange",
    padj_col="padj",
    padj_thresh=0.05,
    lfc_thresh=1.0,
    title=None,
    figsize=(8, 8),
    alpha=0.6,
    point_size=18,
    annotate_top_n=0,
):
    if gene_col is None:
        if df_x.index.name is not None and df_y.index.name is not None:
            x = df_x[[lfc_col, padj_col]].copy().reset_index()
            y = df_y[[lfc_col, padj_col]].copy().reset_index()
            gene_col = df_x.index.name
        else:
            raise ValueError("Provide gene_col or use dataframes indexed by gene IDs.")
    else:
        x = df_x[[gene_col, lfc_col, padj_col]].copy()
        y = df_y[[gene_col, lfc_col, padj_col]].copy()

    x = x.rename(columns={lfc_col: f"{label_x}_lfc", padj_col: f"{label_x}_padj"})
    y = y.rename(columns={lfc_col: f"{label_y}_lfc", padj_col: f"{label_y}_padj"})

    merged = x.merge(y, on=gene_col, how="inner")

    merged["x_sig"] = (
        merged[f"{label_x}_padj"].notna()
        & (merged[f"{label_x}_padj"] < padj_thresh)
        & (merged[f"{label_x}_lfc"].abs() >= lfc_thresh)
    )
    merged["y_sig"] = (
        merged[f"{label_y}_padj"].notna()
        & (merged[f"{label_y}_padj"] < padj_thresh)
        & (merged[f"{label_y}_lfc"].abs() >= lfc_thresh)
    )

    def quadrant(row):
        xlfc = row[f"{label_x}_lfc"]
        ylfc = row[f"{label_y}_lfc"]
        if xlfc >= 0 and ylfc >= 0:
            return "++"
        if xlfc < 0 and ylfc >= 0:
            return "-+"
        if xlfc < 0 and ylfc < 0:
            return "--"
        return "+-"

    merged["quadrant"] = merged.apply(quadrant, axis=1)

    colors = {
        "++": "#d62728",
        "--": "#1f77b4",
        "+-": "#ff7f0e",
        "-+": "#2ca02c",
    }

    fig, ax = plt.subplots(figsize=figsize)

    for q in ["++", "--", "+-", "-+"]:
        sub = merged[merged["quadrant"] == q]
        ax.scatter(
            sub[f"{label_x}_lfc"],
            sub[f"{label_y}_lfc"],
            s=point_size,
            alpha=alpha,
            color=colors[q],
            label=q,
            edgecolors="none",
        )

    ax.axhline(0, color="black", lw=1, alpha=0.7)
    ax.axvline(0, color="black", lw=1, alpha=0.7)

    lim = np.nanmax(
        np.abs(
            merged[[f"{label_x}_lfc", f"{label_y}_lfc"]].to_numpy(dtype=float)
        )
    )
    if np.isfinite(lim):
        ax.set_xlim(-lim * 1.05, lim * 1.05)
        ax.set_ylim(-lim * 1.05, lim * 1.05)

    ax.set_xlabel(f"{label_x} log2FC")
    ax.set_ylabel(f"{label_y} log2FC")
    ax.set_title(title or f"{label_x} vs {label_y}")
    ax.legend(title="Quadrant", loc="best", frameon=False)

    if annotate_top_n > 0:
        top = merged.assign(
            score=merged[f"{label_x}_lfc"].abs() + merged[f"{label_y}_lfc"].abs()
        ).nlargest(annotate_top_n, "score")
        for _, r in top.iterrows():
            ax.text(
                r[f"{label_x}_lfc"],
                r[f"{label_y}_lfc"],
                str(r[gene_col]),
                fontsize=8,
                alpha=0.9,
            )

    plt.tight_layout()
    return fig, ax, merged