## Wrangling Scripts

This file provide a brief overview of the wrangler scripts, and their usage.

### migrations/scripts/migrate_dcp_one_to_dcp_two.py

This script reads the files that were a part of a single DCP 1.0 project, 
transforms them into valid DCP 2.0 metadata objects, and inputs them into 
the DCP 2.0 ledger.  Specify the location to the DCP 1.0 files, and the number
of threads to use for accessing S3.

Example use:
```
    python3 migrations/scripts/migrate_dcp_one_to_dcp_two.py -i <s3dir> --threads <nthreads>
```

Concrete example:
```
    python3 migrations/scripts/migrate_dcp_one_to_dcp_two.py \
       -i  s3://dcp-ledger-bucket-dev/single_cell_transcriptome_analysis_of_human_pancreas_metadata/
       --threads 20
```

