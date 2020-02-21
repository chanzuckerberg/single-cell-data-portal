# Scripts to copy data from DCP 1.0

### migrations/scripts/get_project_metadata.py


This script will download bundles, bundle_manifests and files and store them where the function was called for all of 
the projects in the project_list
```bash
├── hold_data
│   ├── ProjectShortName/
│   │   ├── bundle_manifests/
│   │   ├── bundles/
│   │   ├── errors
│   │   └── files/
```
Example use
```bash
python3 ./dcp_prototype/backend/wrangling/migrations/scripts/get_project_metadata.py -p <project_list> --threads <nthreads>
```
Concrete use
```bash
python3 ./dcp_prototype/backend/wrangling/migrations/scripts/get_project_metadata.py -p ./dcp_prototype/backend/wrangling/migrations/projects.json --threads 20
```
If --thread not passed in, defaults to 1

### migrations/scripts/load_project_metadata.py

This script will load all files (without repeats within a project) and will rename files to {uuid}.{version}.json (with exception of links.json which will be renamed to its checksum)
Example use
```bash
python3 ./dcp_prototype/backend/wrangling/migrations/scripts/load_project_metadata.py -i <input_directory> -o <s3dir> -p <project_list> --threads 20
```
Concrete Use
```bash
python3 ./dcp_prototype/backend/wrangling/migrations/scripts/load_project_metadata.py -i ./hold_data -o https://hca-dcp-one-backup-data.s3.amazonaws.com/dcp-one-copy -p ./dcp_prototype/backend/wrangling/migrations/projects.json -t 20
```
- Can also be used with local directory by passing path with -o arg
- If --thread not passed in, defaults to 1
- The script will search the input_directory `-i` for directories with project names in project_list `-p`

### migrations/scripts/copy_data_files.py
This script will read in the bundle_manifests and compute the blobs (concatenation of hashes) to copy the data files 
from one s3 bucket to another for every project in the input directory 
Example Use
```bash
python3 ./dcp_prototype/backend/wrangling/migrations/scripts/copy_data_files.py -i <input_directory>
```
Concrete Use
```bash
python3 ./dcp_prototype/backend/wrangling/migrations/scripts/copy_data_files.py -i ./hold_data
```
Buckets are currently hard coded:
SOURCE_BUCKET = 'org-hca-dss-prod'
DESTINATION_BUCKET = 'dunitz-prod-copy'


