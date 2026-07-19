from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import yaml


@dataclass
class ProjectConfig:
    fit_dir: str = "./deseq2_fit"
    out_dir: str = "./deseq2_fit/deseq2_results"
    model_file: str = "./deseq2_fit/fitted_dds_anndata.pkl"
    counts_file: str = "./gene-counts/gene_level_counts.csv"
    metadata_file: str = "./sample_metadata.csv"
    log_dir: str = "./logs"

    genotype_col: str = "genotype"
    treatment_col: str = "treatment"
    reference_genotype: str = "WT"
    mutant_genotype: str = "Y537S"
    reference_treatment: str = "shGFP"
    alt_treatment: str = "shPR"

    padj_threshold: float = 0.05
    log2fc_threshold: float = 1.0
    min_total_count: int = 10
    n_cpus: int = 4

    contrasts: dict[str, list[float]] = field(default_factory=dict)

    @classmethod
    def from_yaml(cls, path: str | Path) -> "ProjectConfig":
        path = Path(path)
        with path.open("r") as f:
            data: dict[str, Any] = yaml.safe_load(f) or {}
        return cls(**data)

    def to_dict(self) -> dict[str, Any]:
        return {
            "fit_dir": self.fit_dir,
            "out_dir": self.out_dir,
            "model_file": self.model_file,
            "counts_file": self.counts_file,
            "metadata_file": self.metadata_file,
            "log_dir": self.log_dir,
            "genotype_col": self.genotype_col,
            "treatment_col": self.treatment_col,
            "reference_genotype": self.reference_genotype,
            "mutant_genotype": self.mutant_genotype,
            "reference_treatment": self.reference_treatment,
            "alt_treatment": self.alt_treatment,
            "padj_threshold": self.padj_threshold,
            "log2fc_threshold": self.log2fc_threshold,
            "min_total_count": self.min_total_count,
            "n_cpus": self.n_cpus,
            "contrasts": self.contrasts,
        }
