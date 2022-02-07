import os


def fixture_file_path(filename):
    return os.path.abspath(os.path.join(os.path.dirname(__file__), filename))
