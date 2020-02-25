# Scripts to copy data from DCP 1.0
To download data from dcp 1.0 
- Create a json file containing a list of project shortnames 
- Call`get_project_metadata` passing in the project list
- Load the metadata files into another S3 bucket (or a different directory locally) by calling `load_project_metadata` with the project_list, input_director and output_directory
- Copy the data files between S3 buckets by calling `copy_data_files` 


### migrations/scripts/get_project_metadata.py


This script will download bundles, bundle_manifests and files and store them where the function was called for all of the projects in the project_list
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

This script will load all files in all bundles nested within {input_directory}/{project_name}/bundles/{bundle_id}/ (without repeats within a project) for all projects listed in the passed in project list. Files will be renamed to  {uuid}.{version}.json (with exception of links.json which will be renamed to its checksum)
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
This script will read in the bundle_manifests and compute the blobs (concatenation of hashes) to copy the data files from one s3 bucket to another for every project in the input directory 
Example Use
```bash
python3 ./dcp_prototype/backend/wrangling/migrations/scripts/copy_data_files.py -i <input_directory> -s <source_bucket> -d <destination_bucket>
```
Concrete Use
```bash
python3 ./dcp_prototype/backend/wrangling/migrations/scripts/copy_data_files.py -i ./hold_data -s org-hca-dss-prod -d hca-dcp-one-backup-data
```


