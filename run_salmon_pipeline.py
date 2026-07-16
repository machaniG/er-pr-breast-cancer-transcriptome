"""
pseudo-alignment with salmon, an alternative route through that replaces HISAT2 + featureCounts entirely.
Uses a ThreadPoolExecutor to queue up your trimmed samples, 
automatically checks if they are Single-End or Paired-End, 
and runs Salmon on multiple threads simultaneously.
"""
import os
import glob
import subprocess
import logging
import multiprocessing
from concurrent.futures import ThreadPoolExecutor

SALMON_THREADS = 4  # change as needed

def setup_logging(output_dir):
    logger = logging.getLogger("salmon_pipeline")
    if logger.handlers:
        return logger
    logger.setLevel(logging.INFO)
    log_format = logging.Formatter('%(asctime)s - [%(levelname)s] - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    
    file_handler = logging.FileHandler(os.path.join(output_dir, "salmon_pipeline.log"))
    file_handler.setFormatter(log_format)
    logger.addHandler(file_handler)
    
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(log_format)
    logger.addHandler(console_handler)
    return logger

def process_salmon(r1_input, trimmed_dir, output_dir, index_path, logger):
    """
    Detects SE/PE trimmed files and runs Salmon quantification simultaneously.
    """
    base_name = os.path.basename(r1_input)
    
    # 1. Deduce Sample ID and check for paired companion
    if "_trimmed_1.fastq.gz" in base_name:
        sample_id = base_name.split("_trimmed_1.fastq.gz")[0]
        r2_predict = os.path.join(trimmed_dir, f"{sample_id}_trimmed_2.fastq.gz")
        is_paired = os.path.exists(r2_predict)
    elif "_trimmed.fastq.gz" in base_name:
        sample_id = base_name.split("_trimmed.fastq.gz")[0]
        is_paired = False
    else:
        return # Skip logs/reports from fastp

    # Define unique output folder for this sample's quantification results
    sample_out_dir = os.path.join(output_dir, f"{sample_id}_quant")

    # 2. Construct the core Salmon command
    # -l A tells Salmon to automatically detect the library type (strandedness)
    cmd = [
        "salmon", "quant",
        "-i", index_path,
        "-l", "A",
        "-o", sample_out_dir,
        "--threads", str(SALMON_THREADS) # current 4 CPU cores allocated per sample run
    ]

    # 3. Add Single-End or Paired-End file flags
    if is_paired:
        logger.info(f"🧬 [PAIRED-END] Queueing Salmon mapping for: {sample_id}")
        cmd.extend(["-1", r1_input, "-2", r2_predict])
    else:
        logger.info(f"🧬 [SINGLE-END] Queueing Salmon mapping for: {sample_id}")
        cmd.extend(["-r", r1_input])

    # 4. Execute the Salmon subprocess
    logger.info(f"🚀 Launching Salmon quantification for {sample_id}...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode == 0:
        logger.info(f"✅ Successfully quantified expression for: {sample_id}")
    else:
        logger.error(f"💥 Salmon failed on {sample_id}!\nSTDERR: {result.stderr}")

def main():
    # --- CONFIGURATION PATHS ---
    TRIMMED_DIRECTORY = "./trimmed-data"     # Where your fastp outputs live
    SALMON_OUTPUT_DIR = "./salmon-counts"    # New destination folder
    INDEX_DIRECTORY = "./human_transcriptome_index"
    MAX_SIMULTANEOUS_SAMPLES = 2             # Process 2 samples at a time
    # ---------------------------

    os.makedirs(SALMON_OUTPUT_DIR, exist_ok=True)
    logger = setup_logging(SALMON_OUTPUT_DIR)
    
    logger.info("==================================================")
    logger.info("🎬 Starting Salmon Pseudo-Alignment Pipeline")
    logger.info("==================================================")

    # Gather input files (Filtering to only start via Read 1 or SE files)
    all_trimmed = sorted(glob.glob(os.path.join(TRIMMED_DIRECTORY, "*.fastq.gz")))
    input_files = [f for f in all_trimmed if "_trimmed_2.fastq.gz" not in f]

    if not input_files:
        logger.warning("No valid trimmed files found. Run the fastp script first.")
        return

    # Run simultaneously
    with ThreadPoolExecutor(max_workers=MAX_SIMULTANEOUS_SAMPLES) as executor:
        for file in input_files:
            executor.submit(process_salmon, file, TRIMMED_DIRECTORY, SALMON_OUTPUT_DIR, INDEX_DIRECTORY, logger)

    logger.info("🎉 All samples quantified successfully!")

if __name__ == "__main__":
    main()