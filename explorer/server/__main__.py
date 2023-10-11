import logging

from server.common.utils.utils import import_plugins

# Work around bug https://github.com/pallets/werkzeug/issues/461
if __package__ is None:
    import sys
    from pathlib import Path

    PKG_PATH = Path(__file__).parent
    sys.path.insert(0, str(PKG_PATH.parent))
    import server  # noqa F401

    __package__ = PKG_PATH.name

# Main thing
from .cli.cli import cli  # noqa F402

try:
    import_plugins("server.plugins")
except Exception as e:
    # Make sure to exit in this case, as the server may not be configured as expected.
    logging.critical(f"Error in import_plugins: {str(e)}")
    sys.exit(1)

cli()
