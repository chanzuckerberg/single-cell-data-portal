import os


def check_file(file_path: str) -> bool:
    """
    Check if a file exists, is not a directory, is readable, and has content.

    Args:
        file_path (str): The path to the file to be checked.

    Returns:
        bool: True if the file passes all checks, False otherwise.
    """
    if os.path.isfile(file_path) and os.access(file_path, os.R_OK) and os.path.getsize(file_path) > 0:
        try:
            with open(file_path, "rb") as f:
                content = f.read(1024)
                if content:
                    return True
        except Exception:
            return False
    return False
