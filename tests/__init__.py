import os
import sys

os.environ["CORPORA_HOME"] = os.path.join(os.path.dirname(__file__), "..")  # noqa
pkg_root = os.path.abspath(os.environ["CORPORA_HOME"])  # noqa
sys.path.insert(0, pkg_root)  # noqa
