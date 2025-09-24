import os


def check_file(file_path: str) -> bool:
    """
    Check if a file exists, is not a directory, is readable, and has content.

    Args:
        file_path (str): The path to the file to be checked.

    Returns:
        bool: True if the file passes all checks, False otherwise.
    """
    return os.path.isfile(file_path) and os.path.getsize(file_path) > 1
