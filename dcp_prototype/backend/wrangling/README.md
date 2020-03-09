# Wrangling

Currently the wrangling directory primarily consists of scripts that are required to convert all DCP 1.0 data into a format that allows for its consumption by DCP 2.0. 

## Migration Scripts

All migration scripts, primarily written for the purpose of converting DCP 1.0 data into a format suitable for the DCP 2.0 backend, are hosted in `migrations/scripts`. 

### `migrate_dcp_one_to_dcp_two.py`

This script reads the metadata JSON files that were a part of a single DCP 1.0 project, transforms them into the new artifact schema, and outputs the created JSON file. 

The location of the DCP 1.0 files can be one of the following:
 * an S3 bucket containing metadata files as individual objects
 * an S3 object which contains all DCP 1.0 metadata as a `tar.gz` file 
 * a local `tar.gz` file
 * a local directory containing individual files
 
#### Performance

Reading from a `tar.gz` file is much faster than reading either an S3 bucket with the individual file or reading from a local directory containing each of the metadata files separately. If one of the latter two options is chosen, it is recommended that the script be run with multiple threads. Using a `tar.gz` file to generate the artifact-specified JSON file output does not require multiple threads.

#### Usage

If specified, the `--output-file` is the location where the json is saved. If not specified, the output goes to standard out.

If `--output-file` is specified, and `--append` is specified, then the new project will be combined with the existing output file.  It will add a new project if it is new, or replace an existing project if it is repeated.

When specified, `--threads` allows you to specify the number of threads used to process the metadata files.

Example use:
```
python3 migrations/scripts/migrate_dcp_one_to_dcp_two.py -i <source> --threads <nthreads> --output-file <output> --append
```

Concrete example:
```
# Location in S3 where the DCP projects are found
export BUCKET=s3://hca-dcp-one-backup-data/dcp-one-copy

# Name of the DCP project to process
export PROJECT=Single-cell-transcriptome-analysis-of-human-pancreas

# Run script
python3 migrations/scripts/migrate_dcp_one_to_dcp_two.py -i $BUCKET/$PROJECT.tar.gz  --output-file artifact.json --append
```

### `metadata_aggregator.py`

Given a single DCP 1.0 project's metadata that is contained in an S3 bucket, this script packages it all into a single `tar.gz` file and uploads the created package back up to the S3 bucket. In the process of this transformation, the script will download all of the metadata files locally.

#### Usage

`-i` is the input directory; this should be an S3 bucket with all the metadata files associated with a single project.

`-d` is the input directory; this allows you to specify a local directory from which the metadata files associated with a single project should be packaged into the `tar.gz`.

`--threads` allows you to specify the number of threads used to process the metadata files.

`-o` is the output directory which is the location to which the `tar.gz` package should be uploaded.

Example use:
```
python3 migrations/scripts/metadata_aggregator.py -i <s3dir> -d <localdir> --threads <nthreads> -o <s3output>
```

Running the above command will pull files from `s3dir`, copy them into `localdir` (using `nthreads`), and then will upload the resulting `tar.gz` file back to S3 at the `s3output` location.

Concrete example:
```
# Location in S3 where the DCP projects are found
export BUCKET=s3://hca-dcp-one-backup-data/dcp-one-copy

# Name of the DCP project to process
export PROJECT=Single-cell-transcriptome-analysis-of-human-pancreas

# Create a local directory to which the project should be locally downloaded
mkdir $PROJECT

# Run script
python3 migrations/scripts/metadata_aggregator.py -i $BUCKET/$PROJECT -d $PROJECT --threads 10 -o $BUCKET/$PROJECT.tar.gz
``` 

*** NOTE: This script has already been run on all the 29 DCP 1.0 projects, therefore it may no longer be needed. However if the project data changes, this may potentially need to be rerun. Although if the project metadata needs to be regenerated, then it is simpler to produce the `tar.gz` locally and then upload that file.

### `create_artifact.py`

This script creates an artifact from all the DCP 1.0 projects. By default it chooses all the `.tar.gz` files in `s3://hca-dcp-one-backup-data/dcp-one-copy/`.
That location can be changed with a command line parameter.

#### Usage

`--input_directory` or `-i` allows you to specify the S3 bucket and key containing all of the DCP 1.0 projects.

`--output-file` or `-o` allows you to specify the name of the JSON file that will be outputted, formatted to the Artifact specification.

Example use:
```
python3 migrations/scripts/create_artifact.py [-i <s3 source>] -o <output json file>
```

Concrete example:
```
python3 migrations/scripts/create_artifact.py -o Artifact.json
```
