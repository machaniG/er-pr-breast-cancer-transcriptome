# root/collapse_counts.py
"""
collapse_counts.py — Collapse Salmon transcript-level quantifications
                      into a gene-level count matrix for PyDESeq2.

REFERENCE FORMAT
  This version is built for GENCODE-style transcript FASTA references
  (e.g. gencode.v50.transcripts.fa.gz), where each header is a
  pipe-delimited string:

      >ENST00000456328.2|ENSG00000223972.5|OTTHUMG...|...|DDX11L1-202|DDX11L1|1657|
         field 0 = transcript_id (versioned)
         field 1 = gene_id       (versioned)

  Salmon preserves the full header (up to the first whitespace) as the
  "Name" column in quant.sf. Since GENCODE headers have no internal
  whitespace, the WHOLE pipe-delimited string becomes the Name. We
  exploit this: the gene ID is already embedded in every transcript's
  name, so no external tx2gene download is needed — and no risk of a
  reference-version mismatch, since the mapping is self-consistent by
  construction.

  If you ever switch to a plain Ensembl cDNA FASTA (headers like
  ">ENST00000456328.2 cdna chromosome:..."), this extraction method
  won't apply — you'd go back to an external tx2gene map instead.

WHAT THIS DOES
    1. Locates all sample quant.sf files
    2. Validates the reference ID format on one file (fail fast)
    3. Converts TPM -> length-corrected scaled counts (tximport's
       "scaledTPM" method — removes transcript-length bias)
    4. Extracts each transcript's gene ID and sums counts per gene
    5. Combines all samples into one genes x samples matrix
    6. Saves it to a subfolder

LOGGING
  Every run prints color-coded progress to the console AND writes a
  permanent, timestamped log file to ./logs/.

HOW TO RUN
    cd /path/to/your/project_root
    python3 collapse_counts.py

REQUIREMENTS (one-time install)
    pip install pandas
"""

import glob
import logging
import os
import re
import sys
from dataclasses import dataclass
from datetime import datetime

import pandas as pd


# ==============================================================
# CONFIGURATION
# ==============================================================
@dataclass
class Config:
    salmon_glob: str = "./salmon-counts/*_quant/quant.sf"  # where Salmon results live
    output_dir: str = "./gene-counts"                       # subfolder for results
    output_filename: str = "gene_level_counts.csv"
    log_dir: str = "./logs"                                  # subfolder for run logs
    # Below this %, something is wrong with the ID format assumption.
    min_acceptable_mapping_rate: float = 80.0


CONFIG = Config()


# ==============================================================
# LOGGING SETUP
# ==============================================================
logger = logging.getLogger(__name__)


class ColorFormatter(logging.Formatter):
    COLORS = {
        logging.INFO: "\033[0;32m",
        logging.WARNING: "\033[1;33m",
        logging.ERROR: "\033[0;31m",
        logging.CRITICAL: "\033[1;31m",
    }
    RESET = "\033[0m"

    def format(self, record):
        color = self.COLORS.get(record.levelno, "")
        message = super().format(record)
        return f"{color}{message}{self.RESET}"


def setup_logging(log_dir: str) -> str:
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(
        log_dir, f"collapse_counts_{datetime.now():%Y%m%d_%H%M%S}.log"
    )

    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(ColorFormatter("  [%(levelname)-7s] %(message)s"))

    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s  %(levelname)-7s  %(message)s")
    )

    logger.addHandler(console_handler)
    logger.addHandler(file_handler)
    return log_path


def section(title: str) -> None:
    logger.info("")
    logger.info("─" * 60)
    logger.info(title)


def die(msg: str) -> None:
    logger.error(msg)
    sys.exit(1)


# ==============================================================
# FUNCTION: strip_id_version
# ==============================================================
def strip_id_version(transcript_id: str) -> str:
    return re.sub(r"\.\d+$", "", str(transcript_id))


# ==============================================================
# FUNCTION: find_quant_files
# ==============================================================
def find_quant_files(pattern: str) -> list:
    quant_files = sorted(glob.glob(pattern))

    if not quant_files:
        die(
            f"No quant.sf files found matching: {pattern} | "
            f"Check that you're running this from your project root, "
            f"and that Salmon has already finished (Step 7)."
        )

    logger.info(f"Found {len(quant_files)} sample(s):")
    for f in quant_files:
        sample = os.path.basename(os.path.dirname(f))
        logger.debug(f"  - {sample}")

    return quant_files


# ==============================================================
# FUNCTION: validate_reference_format
# ==============================================================
# PURPOSE : Peek at the FIRST quant.sf file only, before processing
#           all 11 samples, and confirm the Name column actually has
#           the expected GENCODE pipe-delimited shape. This is a
#           fail-fast check — it would catch the 0% mapping
#           issue immediately on one file instead of after a full run.
# ==============================================================
def validate_reference_format(file_path: str) -> None:
    sample_df = pd.read_csv(file_path, sep="\t", nrows=3)
    example_name = sample_df["Name"].iloc[0]
    fields = example_name.split("|")

    logger.info(f"Checking reference ID format using: {os.path.basename(file_path)}")
    logger.info(f"  Example Name value: {example_name}")
    logger.info(f"  Pipe-delimited fields found: {len(fields)}")

    if len(fields) < 2:
        die(
            "This doesn't look like a GENCODE-style pipe-delimited transcript name "
            "(expected 'transcript_id|gene_id|...'). Check what reference FASTA "
            "your Salmon index was actually built from — this script assumes GENCODE."
        )

    example_gene_id = strip_id_version(fields[1])
    logger.info(f"  Extracted gene_id (field 1): {example_gene_id}")


# ==============================================================
# FUNCTION: extract_gene_ids
# ==============================================================
# PURPOSE : Vectorized extraction of gene IDs directly from GENCODE's
#           pipe-delimited transcript names — no external download,
#           no possibility of a reference-version mismatch.
#
#           "ENST00000456328.2|ENSG00000223972.5|...|DDX11L1-202|..."
#                                ^^^^^^^^^^^^^^^^^^  <- field 1, this
# ==============================================================
def extract_gene_ids(names: pd.Series) -> pd.Series:
    fields = names.str.split("|", n=2, expand=True)
    if fields.shape[1] < 2:
        return pd.Series([None] * len(names), index=names.index)
    return fields[1].map(strip_id_version)


# ==============================================================
# FUNCTION: scale_transcript_counts
# ==============================================================
# PURPOSE : Convert TPM -> length-corrected "scaledTPM" counts:
#               raw            = TPM * EffectiveLength
#               scaled_counts  = (raw / sum(raw)) * library_size
# ==============================================================
def scale_transcript_counts(df: pd.DataFrame) -> pd.DataFrame:
    library_size = df["NumReads"].sum()

    raw = df["TPM"] * df["EffectiveLength"]
    raw_sum = raw.sum()

    if raw_sum == 0:
        die("Sample has zero usable signal (TPM * EffectiveLength sums to 0). "
            "Check this sample's Salmon output.")

    df["scaled_counts"] = (raw / raw_sum) * library_size
    return df


# ==============================================================
# FUNCTION: load_sample_gene_counts
# ==============================================================
def load_sample_gene_counts(file_path: str) -> tuple:
    sample_name = os.path.basename(os.path.dirname(file_path)).split("_quant")[0]

    df = pd.read_csv(file_path, sep="\t")
    total_transcripts = len(df)

    df = scale_transcript_counts(df)
    df["gene_id"] = extract_gene_ids(df["Name"])

    mapped_df = df.dropna(subset=["gene_id"])
    mapped_count = len(mapped_df)
    mapping_rate = 100 * mapped_count / total_transcripts if total_transcripts else 0

    gene_counts = mapped_df.groupby("gene_id")["scaled_counts"].sum()

    diagnostics = {
        "sample": sample_name,
        "total_transcripts": total_transcripts,
        "mapped_transcripts": mapped_count,
        "mapping_rate": mapping_rate,
        "genes_recovered": gene_counts.shape[0],
    }

    logger.debug(f"{sample_name}: {mapping_rate:.1f}% mapped, {gene_counts.shape[0]:,} genes")
    return sample_name, gene_counts, diagnostics


# ==============================================================
# FUNCTION: build_count_matrix
# ==============================================================
def build_count_matrix(quant_files: list) -> pd.DataFrame:
    all_counts = {}
    all_diagnostics = []

    for file_path in quant_files:
        sample_name, gene_counts, diag = load_sample_gene_counts(file_path)
        all_counts[sample_name] = gene_counts
        all_diagnostics.append(diag)

    logger.info(f"{'Sample':<20}{'Transcripts':>13}{'Mapped':>10}{'Mapping %':>12}{'Genes':>9}")
    logger.info("-" * 64)
    for d in all_diagnostics:
        row = (
            f"{d['sample']:<20}{d['total_transcripts']:>13,}"
            f"{d['mapped_transcripts']:>10,}{d['mapping_rate']:>11.1f}%"
            f"{d['genes_recovered']:>9,}"
        )
        if d["mapping_rate"] < CONFIG.min_acceptable_mapping_rate:
            logger.warning(f"{row}  <-- LOW MAPPING RATE")
        else:
            logger.info(row)

    low_samples = [d for d in all_diagnostics if d["mapping_rate"] < CONFIG.min_acceptable_mapping_rate]
    if low_samples:
        logger.warning(
            f"{len(low_samples)} sample(s) below {CONFIG.min_acceptable_mapping_rate}% mapping rate. "
            f"Re-run validate_reference_format() and inspect quant.sf 'Name' values directly."
        )
    else:
        logger.info("All samples above minimum mapping rate threshold.")

    final_matrix = pd.DataFrame(all_counts).fillna(0)
    final_matrix = final_matrix.round().astype(int)
    return final_matrix


# ==============================================================
# FUNCTION: save_matrix
# ==============================================================
def save_matrix(matrix: pd.DataFrame, output_dir: str, filename: str) -> str:
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, filename)
    matrix.to_csv(out_path)
    return out_path


# ==============================================================
# MAIN
# ==============================================================
def main():
    log_path = setup_logging(CONFIG.log_dir)

    logger.info("=" * 60)
    logger.info("Salmon Transcript Counts -> Gene-Level Counts (GENCODE)")
    logger.info(f"Log file: {log_path}")
    logger.info("=" * 60)

    try:
        section("1/4  Locating Salmon quant.sf files")
        quant_files = find_quant_files(CONFIG.salmon_glob)

        section("2/4  Validating reference ID format")
        validate_reference_format(quant_files[0])

        section("3/4  Scaling and collapsing to gene level")
        final_matrix = build_count_matrix(quant_files)

        section("4/4  Saving results")
        out_path = save_matrix(final_matrix, CONFIG.output_dir, CONFIG.output_filename)
        logger.info(f"Matrix shape: {final_matrix.shape[0]:,} genes x {final_matrix.shape[1]} samples")
        logger.info(f"Saved to: {out_path}")

    except SystemExit:
        raise
    except Exception:
        logger.exception("Unexpected error — see traceback above and in the log file.")
        sys.exit(1)

    logger.info("")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()

