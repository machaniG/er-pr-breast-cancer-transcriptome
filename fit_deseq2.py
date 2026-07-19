#!/usr/bin/env python3
from __future__ import annotations
import logging

import json
import pickle as pkl
from dataclasses import dataclass, asdict
from pathlib import Path
import sys
import os
import pandas as pd
from pydeseq2.dds import DeseqDataSet

try:
    from pydeseq2.ds import DeseqStats
except Exception:
    DeseqStats = None

logging.basicConfig()
logger = logging.getLogger(__name__)

@dataclass
class Config:
    counts_file: str = "./gene-counts/gene_level_counts.csv"
    metadata_file: str = "./sample_metadata.csv"
    output_dir: str = "./deseq2_fit"
    reference_genotype: str = "WT"
    mutant_genotype: str = "mutant"
    reference_treatment: str = "shGFP"
    alt_treatment: str = "shPR"
    genotype_col: str = "genotype"
    treatment_col: str = "treatment"
    padj_threshold: float = 0.05
    log2fc_threshold: float = 1.0
    min_total_count: int = 10
    n_cpus: int = 4
    model_file: str = "./deseq2_fit/deseq_model.pkl"



def load_counts(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, index_col="gene_id")
    #filter low reads
    df = df[df.sum(axis=1) >= 10]
    # transpose
    df = df.T

    return df


def load_metadata1(path: str, sample_names: list[str]) -> pd.DataFrame:
    md = pd.read_csv(path, index_col='sample_id')
    md = md.loc[sample_names].copy()
    md["genotype"] = pd.Categorical(md["genotype"], categories=["WT", "mutant"], ordered=True)
    md["treatment"] = pd.Categorical(md["treatment"], categories=["shGFP", "shPR"], ordered=True)
    return md

def load_metadata(path: str) -> pd.DataFrame:
    md = pd.read_csv(path, index_col = "sample_id")

    return md

def save_deseq_model(dds: DeseqDataSet, config: Config) -> Path:
    model_path = Path(config.model_file)
    model_path.parent.mkdir(parents=True, exist_ok=True)
    with open(model_path, "wb") as fh:
        pkl.dump(dds, fh, protocol=pkl.HIGHEST_PROTOCOL)
    logger.info("Saved fitted DESeq2 model object to %s", model_path)
    return model_path

def main():
    cfg = Config()
    outdir = Path(cfg.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    counts = load_counts(cfg.counts_file)
    md = load_metadata(cfg.metadata_file)

    shared = counts.index.intersection(md.index)
    counts_df = counts.loc[shared].sort_index()
    metadata_df = md.loc[shared].sort_index()

    assert counts_df.index.equals(metadata_df.index), "Sample IDs do not match!"
    # sort both dataframes so the sample IDs match perfectly row for row
    #counts_df = counts.sort_index()
    #metadata_df = md.sort_index()
    # Double-check that they match exactly
    #assert (counts_df.index == metadata_df.index).all(), "Sample IDs do not match!"

    
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata_df,
        design="~ genotype + treatment + genotype:treatment",
        n_cpus=cfg.n_cpus,
    )
    dds.deseq2()

    # 3. Convert to a safe format and save to disk
    os.makedirs("deseq2_fit", exist_ok=True)
    with open("deseq2_fit/fitted_dds_anndata.pkl", "wb") as f:
        pkl.dump(dds.to_picklable_anndata(), f)
    
    logger.info("Model fitting complete and object saved successfully!")

    with open(outdir / "metadata_used.csv", "w") as f:
        metadata_df.to_csv(f)
    logger.info("Metadata saved successfully!")

    # Convert DeseqStats results_df to CSV
    #stat_res.results_df.to_csv(os.path.join(output_path, "diff_expression_results.csv"))

    #Print the raw column names
    print("\n--- Model Coefficients Found ---")
    coef_names = list(dds.varm["LFC"].columns)
    print(coef_names)

    

if __name__ == "__main__":
    main()