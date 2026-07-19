"""
custom logging set up 
use:
from config import ProjectConfig
from logging_utils import setup_logging

cfg = ProjectConfig.from_yaml("config.yml")
logger = setup_logging(Path(cfg.out_dir) / "pipeline.log")
change the last argument for each script to get a separate log file for each run
e.g.,  
        logger = setup_logging(Path(cfg.out_dir) / "fit_deseq2.log") 
        logger = setup_logging(Path(cfg.out_dir) / "alignment.log")
"""

from __future__ import annotations

import logging
from pathlib import Path


def setup_logging(log_file: str | Path, logger_name: str = "rnaseq_pipeline") -> logging.Logger:
    log_file = Path(log_file)
    log_file.parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    logger.propagate = False

    if logger.handlers:
        return logger

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    return logger

