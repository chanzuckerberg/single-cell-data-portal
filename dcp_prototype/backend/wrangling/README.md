## Wrangling Scripts

This file provide a brief overview of the wrangler scripts, and their usage.

### migrations/scripts/migrate_dcp_one_to_dcp_two.py

This script reads the files that were a part of a single DCP 1.0 project, 
transforms them into the new artifact schema, and outputs a json file. 

The location can be one of the following:
 * an s3 bucket containing metadata files and individual object
 * an s3 object which contains all dcp1 metadata as a tar.gz file 
 * a local tar.gz file
 * a local directory containing individual files
 
It is much much faster to use a tar.gz file.
When using a .tar.gz file it is no longer necessary to run with multiple threads.

If specified, the --output-file is the location where the json is saved.
If not specified, the output goes to standard out.

If --output-file is specified, and --append is specified, then the new project
will be combined with the existing output file.  It will add a new project if it is new,
or replace an existing project if it is repeated.

Example use:
```
    python3 migrations/scripts/migrate_dcp_one_to_dcp_two.py -i <source> --threads <nthreads> --output-file <output> --append
```

Concrete example:
```
    # location in S3 where the DCP projects are found
    export BUCKET=s3://hca-dcp-one-backup-data/dcp-one-copy
    # name of the DCP project to process
    export PROJECT=Single-cell-transcriptome-analysis-of-human-pancreas
    python3 migrations/scripts/migrate_dcp_one_to_dcp_two.py -i $BUCKET/$PROJECT.tar.gz  --output-file artifact.json --append
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
    # location in S3 where the DCP projects are found
    export BUCKET=s3://hca-dcp-one-backup-data/dcp-one-copy
    # name of the DCP project to process
    export PROJECT=Single-cell-transcriptome-analysis-of-human-pancreas
    mkdir $PROJECT
    python3 migrations/scripts/metadata_aggregator.py -i $BUCKET/$PROJECT -d $PROJECT --threads 10 -o $BUCKET/$PROJECT.tar.gz 
``` 

This script has already been run on all the 29 DCP 1.0 projects, therefore it may no longer be needed.
However if the project data changes, then potentially this may need to be rerun.  
Although if the project metadata needs to be regenerated, then it is simpler to produce the tar.gz 
locally and then upload that file.

### migrations/scripts/create_artifact.py

Create an artifact from all the dcp1 projects.
By default it chooses all the .tar.gz files in s3://hca-dcp-one-backup-data/dcp-one-copy/ .
That location can be changes with a command line parameter.

Example use:
```
    python3 migrations/scripts/create_artifact.py  [-i <s3 source>] -o <output json file>
```

Concrete example:
```
    python3 migrations/scripts/create_artifact.py  -o Artifact.json
```
