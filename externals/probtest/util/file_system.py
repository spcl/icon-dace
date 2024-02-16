from pathlib import Path

from util.log_handler import logger


def file_names_from_pattern(dir_name, file_pattern):
    """
    Search for all file patching file_pattern in directory dir_name

    Parameters
    ----------
    dir_name : str
    file_pattern : str
        file_pattern may contain simple shell-style wildcards such as "*".
        file_pattern may contain "/" to specify subdirectories.

    Returns
    -------
    file_names : list of str
        A list of matched file names that are relative to dir_name.
    err : int
        non-zero if error occurred
    """
    err = 0
    matched_pattern = Path(dir_name).glob(file_pattern)
    # use .relative_to to allow file_pattern to contain subdirectories.
    file_names = [p.relative_to(dir_name) for p in matched_pattern]
    file_names = [str(f) for f in file_names]
    if len(file_names) < 1:
        logger.warning(
            "no files found in '{}' for pattern '{}'".format(dir_name, file_pattern)
        )
        logger.debug("the directory contains the following files:")
        for f in Path(dir_name).glob("*"):
            logger.debug(f.name)
        err = 1

    return file_names, err
