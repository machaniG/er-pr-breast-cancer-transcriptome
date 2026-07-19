#!/usr/bin/env python3
from __future__ import annotations

import json
import logging
import os
import yaml
from dataclasses import dataclass, asdict
from pathlib import Path
import pickle as pkl
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from g_config import ProjectConfig
from logging_utils import setup_logging

cfg = ProjectConfig.from_yaml("config.yml")
logger = setup_logging(Path(cfg.log_dir) / "extract_contrasts.log", __name__)


def annotate_results(df: pd.DataFrame, cfg: ProjectConfig) -> pd.DataFrame:
    out = df.copy()
    out["significant"] = (
        out["padj"].notna()
        & (out["padj"] < cfg.padj_threshold)
        & (out["log2FoldChange"].abs() >= cfg.log2fc_threshold)
    )
    out["direction"] = "not_significant"
    out.loc[out["significant"] & (out["log2FoldChange"] > 0), "direction"] = "up"
    out.loc[out["significant"] & (out["log2FoldChange"] < 0), "direction"] = "down"
    return out


def save_tables(df: pd.DataFrame, prefix: str, outdir: Path) -> dict[str, Path]:
    df = df.sort_values("padj", na_position="last")
    paths = {
        "all": outdir / f"{prefix}.csv",
        "significant": outdir / f"{prefix}.significant.csv",
        "up": outdir / f"{prefix}.up.csv",
        "down": outdir / f"{prefix}.down.csv",
    }

    df.to_csv(paths["all"])
    df[df["significant"]].to_csv(paths["significant"])
    df[df["direction"] == "up"].to_csv(paths["up"])
    df[df["direction"] == "down"].to_csv(paths["down"])
    return paths


def summarize_result(df: pd.DataFrame) -> dict[str, int]:
    return {
        "tested": int(df["padj"].notna().sum()),
        "significant": int(df["significant"].sum()),
        "up": int((df["direction"] == "up").sum()),
        "down": int((df["direction"] == "down").sum()),
    }


def main():
    #cfg = ProjectConfig()
    fit_dir = Path(cfg.fit_dir)
    outdir = Path(cfg.out_dir)
    #logger = setup_logging(log_dir)

    logger.info("Loading saved model from %s", fit_dir / "fitted_dds_anndata.pkl")
    with open(fit_dir / "fitted_dds_anndata.pkl", "rb") as f:
        loaded_adata = pkl.load(f)

    dds = DeseqDataSet(
        adata=loaded_adata,
        design="~ genotype + treatment + genotype:treatment",
    )

    # print model coeffs    
    try:
        coef_names = list(dds.varm["LFC"].columns)
    except Exception:
        coef_names = []
    logger.info("Model coefficients: %s", coef_names)
    

    """
    with open("config.yml", "r") as f:
        config = yaml.safe_load(f)
   
   
    contrasts_dict = config["contrasts"]
    """

    contrasts_dict = cfg.contrasts

    saved = {}
    for name, vector_list in contrasts_dict.items():
        logger.info("Processing contrast: %s", name)

        contrast = np.array(vector_list, dtype=float)

        stat = DeseqStats(
            dds,
            contrast=contrast,
            alpha=cfg.padj_threshold,
            independent_filter=True,
            cooks_filter=True,
            quiet=True,
        )

        stat.summary()

        res_df = stat.results_df.copy()
        res_df.index.name = "gene_id"
        res_df = annotate_results(res_df, cfg)

        paths = save_tables(res_df, name, outdir)
        saved[name] = {k: str(v) for k, v in paths.items()}

        s = summarize_result(res_df)
        logger.info(
            "%s | tested=%d | significant=%d | up=%d | down=%d",
            name, s["tested"], s["significant"], s["up"], s["down"]
        )


    with open(outdir / "config.json", "w") as f:
        json.dump(asdict(cfg), f, indent=2)

    with open(outdir / "saved_outputs.json", "w") as f:
        json.dump(saved, f, indent=2)

    logger.info("Saved all contrast results to %s", outdir)


if __name__ == "__main__":
    main()

