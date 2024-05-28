import os

MARKER_SCORE_THRESHOLD = 0.5

# If PIPELINE_NUM_CPUS is not set, use 12 CPUs by default
# If the number of CPUs on the machine is less than 12, use the number of CPUs on the machine
# TODO: Tune this number. But note - speed is not that important here and we've run into issues
# where the number of CPUs that was previously set (24) became too high and resulted in OOM isues
# as the data corpus grew.
PIPELINE_NUM_CPUS = min(os.cpu_count(), os.getenv("PIPELINE_NUM_CPUS", 12))
