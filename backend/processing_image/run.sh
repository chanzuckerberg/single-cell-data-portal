#!/bin/bash
set -euo pipefail

ROOT_DIR=/data
ls $ROOT_DIR
SOURCE_ANNDATA=( $ROOT_DIR/*.h5ad )

if [ ! -f "$SOURCE_ANNDATA" ]; then
	echo "Source AnnData file not found!"
	exit 1
fi

VALIDATION_MESSAGE="$(cellxgene schema validate "$SOURCE_ANNDATA")"
VALIDATION_RC=$?

if [ $VALIDATION_RC -ne 0 ]; then
	echo "Validation failed!"
	echo "$VALIDATION_MESSAGE"
	exit $VALIDATION_RC
fi

echo "Extracting metadata:"
python3 extract_metadata.py "$SOURCE_ANNDATA"

Rscript make_seurat.R "$SOURCE_ANNDATA"
CONVERTED_SEURAT="${SOURCE_ANNDATA%.h5ad}".rds
if [ ! -f "$CONVERTED_SEURAT" ]; then
	echo "Seurat conversion failed!"
fi


python3 make_loom.py "$SOURCE_ANNDATA"
CONVERTED_LOOM="${SOURCE_ANNDATA%.h5ad}".loom
if [ ! -f "$CONVERTED_LOOM" ]; then
	echo "Loom conversion failed!"
fi

cellxgene convert -o "${SOURCE_ANNDATA%.h5ad}".cxg -s 10.0 "$SOURCE_ANNDATA"

chmod 666 "${SOURCE_ANNDATA%.h5ad}".rds
chmod 666 "${SOURCE_ANNDATA%.h5ad}".loom
chmod -R 777 "${SOURCE_ANNDATA%.h5ad}".cxg/
