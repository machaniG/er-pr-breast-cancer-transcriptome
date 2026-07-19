from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def prep_volcano_df(csv_path, log2fc_col="log2FoldChange", padj_col="padj"):
    df = pd.read_csv(csv_path, index_col=0).copy()
    df[padj_col] = df[padj_col].fillna(1.0)
    df[log2fc_col] = df[log2fc_col].fillna(0.0)
    df["minus_log10_padj"] = -np.log10(df[padj_col].clip(lower=1e-300))
    return df


def classify_volcano(df, log2fc_col="log2FoldChange", padj_col="padj",
                     padj_threshold=0.05, log2fc_threshold=1.0):
    sig = (df[padj_col] < padj_threshold) & (df[log2fc_col].abs() >= log2fc_threshold)
    up = sig & (df[log2fc_col] > 0)
    down = sig & (df[log2fc_col] < 0)
    ns = ~sig
    return ns, up, down


def plot_two_panel_volcano(
    left_csv,
    right_csv,
    out_png,
    left_title="shPR knockdown effect in WT",
    right_title="shPR knockdown effect in Y537S mutant",
    log2fc_col="log2FoldChange",
    padj_col="padj",
    padj_threshold=0.05,
    log2fc_threshold=1.0,
    figsize=(14, 6),
    point_size=12,
):
    left = prep_volcano_df(left_csv, log2fc_col=log2fc_col, padj_col=padj_col)
    right = prep_volcano_df(right_csv, log2fc_col=log2fc_col, padj_col=padj_col)

    l_ns, l_up, l_down = classify_volcano(
        left, log2fc_col=log2fc_col, padj_col=padj_col,
        padj_threshold=padj_threshold, log2fc_threshold=log2fc_threshold
    )
    r_ns, r_up, r_down = classify_volcano(
        right, log2fc_col=log2fc_col, padj_col=padj_col,
        padj_threshold=padj_threshold, log2fc_threshold=log2fc_threshold
    )

    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)

    def draw_panel(ax, df, ns, up, down, title):
        ax.scatter(df.loc[ns, log2fc_col], df.loc[ns, "minus_log10_padj"],
                   s=point_size, c="#B8B8B8", alpha=0.6, linewidths=0, label="Not significant")
        ax.scatter(df.loc[down, log2fc_col], df.loc[down, "minus_log10_padj"],
                   s=point_size, c="#4C78A8", alpha=0.8, linewidths=0, label="Down")
        ax.scatter(df.loc[up, log2fc_col], df.loc[up, "minus_log10_padj"],
                   s=point_size, c="#D62728", alpha=0.8, linewidths=0, label="Up")

        ax.axvline(-log2fc_threshold, color="black", linestyle="--", linewidth=1)
        ax.axvline(log2fc_threshold, color="black", linestyle="--", linewidth=1)
        ax.axhline(-np.log10(padj_threshold), color="black", linestyle="--", linewidth=1)

        ax.set_title(title)
        ax.set_xlabel("log2 fold change")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    draw_panel(axes[0], left, l_ns, l_up, l_down, left_title)
    draw_panel(axes[1], right, r_ns, r_up, r_down, right_title)

    axes[0].set_ylabel("-log10 adjusted p-value")
    axes[1].legend(frameon=False, loc="upper right")

    fig.suptitle("Synergy vs. Baseline: shPR Knockdown Effects", y=1.02)
    fig.tight_layout()

    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return out_png