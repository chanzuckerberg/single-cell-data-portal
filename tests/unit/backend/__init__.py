from os import popen

PROJECT_ROOT = popen("git rev-parse --show-toplevel").read().strip()
