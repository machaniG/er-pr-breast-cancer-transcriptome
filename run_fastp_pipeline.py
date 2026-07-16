"""
run_fastp_pipeline.py
A script for adapter trimming and filtering. It strips illumina adapters,
removes low-quality bases and reads that are too short
"""

import os
import glob
import subprocess
import logging
from concurrent.futures import ThreadPoolExecutor

def setup_logging(output_dir):
    """
    Configures logging to output to both the terminal (StreamHandler)
    and a persistent log file (FileHandler) inside the output directory.
    """
    log_file = os.path.join(output_dir, "fastp_pipeline.log")
    
    # Create a master logger
    logger = logging.getLogger("fastp_pipeline")
    logger.setLevel(logging.INFO)
    
    # Prevent duplicate handlers if script is re-run in an interactive environment
    if logger.handlers:
        return logger

    # Define a clean, professional log format with timestamps
    log_format = logging.Formatter('%(asctime)s - [%(levelname)s] - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    # 1. File Handler (Writes to trimmed_data/fastp_pipeline.log)
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(log_format)
    logger.addHandler(file_handler)

    # 2. Console Handler (Writes to your VS Code Terminal)
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(log_format)
    logger.addHandler(console_handler)

    return logger

def process_sample(r1_input, output_dir, logger):
    """
    Dynamically detects if a file is Single-End or Paired-End,
    and runs fastp with the appropriate optimized parameters.
    """
    base_name = os.path.basename(r1_input)
    
    # 1. Determine Sample ID and predict possible R2 names
    if "_1.fastq" in base_name:
        r2_predict = r1_input.replace("_1.fastq", "_2.fastq")
        sample_id = base_name.split("_1.fastq")[0]
    elif "_R1_" in base_name:
        r2_predict = r1_input.replace("_R1_", "_R2_")
        sample_id = base_name.split("_R1_")[0]
    else:
        sample_id = base_name.split(".fastq")[0]
        r2_predict = None

    # 2. Check if it's Actively Paired-End or Single-End
    is_paired = r2_predict is not None and os.path.exists(r2_predict)

    # 3. Build Base Command
    cmd = [
        "fastp",
        "-i", r1_input,
        "--qualified_quality_phred", "20",
        "--unqualified_percent_limit", "40",
        "--length_required", "36",
        "--thread", "4"
    ]

    # 4. Append specific layout parameters
    if is_paired:
        logger.info(f"🧬 [PAIRED-END DETECTED] Initializing sample: {sample_id}")
        r1_output = os.path.join(output_dir, f"{sample_id}_trimmed_1.fastq.gz")
        r2_output = os.path.join(output_dir, f"{sample_id}_trimmed_2.fastq.gz")
        
        cmd.extend([
            "-I", r2_predict,
            "-o", r1_output,
            "-O", r2_output,
            "--detect_adapter_for_pe",
            "--correction"
        ])
    else:
        logger.info(f"🧬 [SINGLE-END DETECTED] Initializing sample: {sample_id}")
        r1_output = os.path.join(output_dir, f"{sample_id}_trimmed.fastq.gz")
        
        cmd.extend([
            "-o", r1_output
        ])

    # 5. Add standard report paths
    html_report = os.path.join(output_dir, f"{sample_id}_fastp.html")
    json_report = os.path.join(output_dir, f"{sample_id}_fastp.json")
    cmd.extend(["--html", html_report, "--json", json_report])

    # 6. Execute command
    logger.info(f"🚀 Launching fastp subprocess for {sample_id}...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode == 0:
        logger.info(f"✅ Successfully processed and trimmed: {sample_id}")
    else:
        logger.error(f"💥 Failed processing {sample_id}!\nSTDOUT: {result.stdout}\nSTDERR: {result.stderr}")

def main():
    # --- CONFIGURATION PATHS ---
    INPUT_DIRECTORY = "./raw-data" 
    OUTPUT_DIRECTORY = "./trimmed-data"
    MAX_SIMULTANEOUS_SAMPLES = 2 
    # ---------------------------

    os.makedirs(OUTPUT_DIRECTORY, exist_ok=True)

    # Initialize Logger
    logger = setup_logging(OUTPUT_DIRECTORY)
    logger.info("==================================================")
    logger.info("🎬 Starting Bioinformatics Fastp Trimming Pipeline")
    logger.info("==================================================")

    # Find input files
    all_fastq = sorted(glob.glob(os.path.join(INPUT_DIRECTORY, "*.fastq*")))
    input_files = []
    for f in all_fastq:
        if "_2.fastq" in f or "_R2_" in f:
            continue
        input_files.append(f)

    if not input_files:
        logger.warning(f"No valid FASTQ files found in input directory: '{INPUT_DIRECTORY}'")
        return

    logger.info(f"Found {len(input_files)} distinct biological samples queueing for processing.")

    # Execute simultaneous processing
    with ThreadPoolExecutor(max_workers=MAX_SIMULTANEOUS_SAMPLES) as executor:
        for file in input_files:
            executor.submit(process_sample, file, OUTPUT_DIRECTORY, logger)

    logger.info("🎉 All samples processed. Pipeline execution complete.")

if __name__ == "__main__":
    main()

    # run it: python run_fastp_pipeline.py