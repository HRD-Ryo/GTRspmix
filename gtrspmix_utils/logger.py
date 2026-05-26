
import logging

def setup_logger(logfile):
    logger = logging.getLogger('gtrspmix')
    logger.setLevel(logging.INFO)

    console_handler = logging.StreamHandler()
    file_handler = logging.FileHandler(str(logfile))

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y/%m/%d %H:%M:%S')
    console_handler.setFormatter(formatter)
    file_handler.setFormatter(formatter)

    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

    return logger

