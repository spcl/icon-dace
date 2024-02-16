import logging
import sys

mpl_logger = logging.getLogger("matplotlib")
mpl_logger.setLevel(logging.WARNING)

logger = logging.getLogger()


def initialize_logger(log_level="DEBUG", log_file="probtest.log"):
    formatter = logging.Formatter("%(message)s")

    stream_handler = logging.StreamHandler(stream=sys.stdout)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.setLevel(log_level)

    logger.info("initialized logger with level {}".format(log_level))
