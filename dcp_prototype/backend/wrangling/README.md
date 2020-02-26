## Wrangling Scripts

This file provide a brief overview of the wrangler scripts, and their usage.

### migrations/scripts/migrate_dcp_one_to_dcp_two.py

This script reads the files that were a part of a single DCP 1.0 project, 
transforms them into valid DCP 2.0 metadata objects, and inputs them into 
the DCP 2.0 ledger.  Specify the location to the DCP 1.0 files, and the number
of threads to use for accessing S3.

The location can be one of the following:
 * an s3 bucket containing metadata files and individual object
 * an s3 object which contains all dcp1 metadata as a tar.gz file 
 * a local tar.gz file
 * a local directory containing individual files
 
It is much much faster to use a tar.gz file.
When using a .tar.gz file it is no longer necessary to run with multiple threads.

Example use:
```
    python3 migrations/scripts/migrate_dcp_one_to_dcp_two.py -i <source> --threads <nthreads>
```

Concrete example:
```
    export BUCKET=s3://hca-dcp-one-backup-data/dcp-one-copy
    ecport PROJECT=Single-cell-transcriptome-analysis-of-human-pancreas
    python3 migrations/scripts/migrate_dcp_one_to_dcp_two.py -i $BUCKET/$PROJECT.tar.gz 
```

### migrations/scripts/metadata_aggregator.py

Script to create a tar.gz file from the metadata contained in an S3 bucket.
It first downloads all the files locally, creates the .tar.gz file, then uploads them back to s3.

Example use:
```
    python3 migrations/scripts/metadata_aggregator.py -i <s3dir> -d <localdir> --threads <nthreads> -o <s3output>
```

This will pull files from s3dir, copy them into localdir (using nthreads), and then uploads the resuling tar.gz file back
to s3 at the <s3output> location.

Concrete example:
```
    export BUCKET=s3://hca-dcp-one-backup-data/dcp-one-copy
    ecport PROJECT=Single-cell-transcriptome-analysis-of-human-pancreas
    python3 migrations/scripts/metadata_aggregator.py -i $BUCKET/$PROJECT -d $PROJECT --threads 10 -o $BUCKET/$PROJECT.tar.gz 
``` 

This script has already been run on all the 29 DCP 1.0 projects, therefore it may no longer be needed.
However if the project data changes, then potentially this may need to be rerun.  
Although if the project metadata needs to be regenerated, then it is simpler to produce the tar.gz 
locally and then upload that file.
