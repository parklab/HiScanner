import logging
import sys
from pathlib import Path
from typing import Optional, Union
import colorlog

class CustomFormatter(colorlog.ColoredFormatter):
    """Custom formatter with colors for different log levels"""
    def __init__(self):
        super().__init__(
            "%(log_color)s%(levelname)-8s%(reset)s %(message)s",
            log_colors={
                'DEBUG':    'cyan',
                'INFO':     'black',
                'WARNING': 'yellow',
                'ERROR':   'red',
                'CRITICAL': 'red,bg_white',
            },
            secondary_log_colors={},
            style='%'
        )

def setup_logging(
    log_file: Optional[Union[str, Path]] = None,
    log_level: Union[str, int] = logging.INFO,
    logger_name: str = "hiscanner"
) -> logging.Logger:
    """
    Setup logging configuration for HiScanner.

    Parameters
    ----------
    log_file : str or Path, optional
        Path to log file. If None, logs will only go to console.
    log_level : str or int
        Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    logger_name : str
        Name of the logger instance

    Returns
    -------
    logging.Logger
        Configured logger instance
    """
    # Convert string log level to numeric value if needed
    if isinstance(log_level, str):
        log_level = getattr(logging, log_level.upper())

    # Create logger
    logger = logging.getLogger(logger_name)
    logger.setLevel(log_level)

    # Remove existing handlers
    logger.handlers = []

    # Create console handler with colored output
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(CustomFormatter())
    logger.addHandler(console_handler)

    # Create file handler if log file is specified
    if log_file is not None:
        log_path = Path(log_file)
        
        # Create directory if it doesn't exist
        log_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Create file handler with detailed formatting
        file_handler = logging.FileHandler(log_path)
        file_formatter = logging.Formatter(
            fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)

    return logger

# Create a default logger instance
logger = setup_logging()

def get_logger(name: str = "hiscanner") -> logging.Logger:
    """
    Get a logger instance.

    Parameters
    ----------
    name : str
        Name of the logger instance

    Returns
    -------
    logging.Logger
        Logger instance
    """
    return logging.getLogger(name)

if __name__ == "__main__":
    # Example usage and testing
    logger = setup_logging(log_level="DEBUG")
    logger.debug("This is a debug message")
    logger.info("This is an info message")
    logger.warning("This is a warning message")
    logger.error("This is an error message")
    logger.critical("This is a critical message")