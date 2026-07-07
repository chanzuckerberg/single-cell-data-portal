# Backend Common Dependencies Used by CellGuide Pipeline

This document maps out which code from `backend/common` is being used by `backend/cellguide/pipeline`.

## Summary

The CellGuide pipeline uses code from the following `backend/common` modules:
- `backend.common.census_cube` - Census cube data and utilities
- `backend.common.marker_genes` - Marker gene computation and metadata
- `backend.common.utils` - Utility functions (dataclass, cloudfront, result_notification)
- `backend.common.doi` - DOI cleaning utilities
- `backend.common.providers` - External service providers (Crossref)

---

## Detailed Mapping

### 1. `backend.common.census_cube`

#### 1.1 `backend.common.census_cube.data.snapshot`
**Used in:**
- `backend/cellguide/pipeline/metadata/__init__.py`
- `backend/cellguide/pipeline/source_collections/__init__.py`
- `backend/cellguide/pipeline/ontology_tree/__init__.py`
- `backend/cellguide/pipeline/gpt_descriptions/__init__.py`
- `backend/cellguide/pipeline/computational_marker_genes/__init__.py`
- `backend/cellguide/pipeline/canonical_marker_genes/__init__.py`

**Usage:**
- `sn.load_snapshot()` - Loads the WMG census cube snapshot with a specific schema version
- Used to access cell counts, primary filter dimensions, and other census data

#### 1.2 `backend.common.census_cube.data.ontology_labels`
**Used in:**
- `backend/cellguide/pipeline/ontology_tree/tree_builder.py`
- `backend/cellguide/pipeline/source_collections/source_collections_generator.py`
- `backend/cellguide/pipeline/metadata/metadata_generator.py`
- `backend/cellguide/pipeline/canonical_marker_genes/canonical_markers.py`

**Functions used:**
- `ontology_term_label(term_id)` - Gets the label/name for an ontology term
- `ontology_term_description(term_id)` - Gets the description for an ontology term
- `ontology_term_synonyms(term_id)` - Gets synonyms for an ontology term
- `is_ontology_term_deprecated(term_id)` - Checks if an ontology term is deprecated

#### 1.3 `backend.common.census_cube.utils`
**Used in:**
- `backend/cellguide/pipeline/ontology_tree/tree_builder.py`
- `backend/cellguide/pipeline/source_collections/source_collections_generator.py`
- `backend/cellguide/pipeline/metadata/__init__.py`
- `backend/cellguide/pipeline/canonical_marker_genes/__init__.py`
- `backend/cellguide/pipeline/canonical_marker_genes/canonical_markers.py`

**Functions used:**
- `ancestors(term_id)` - Gets all ancestor terms for an ontology term
- `descendants(term_id)` - Gets all descendant terms for an ontology term
- `ontology_parser` - Parser for ontology data
- `rollup_across_cell_type_descendants(df)` - Rolls up cell counts across descendant cell types
- `to_dict(keys, values)` - Converts two arrays to a dictionary
- `get_all_cell_type_ids_in_corpus(snapshot)` - Gets all cell type IDs present in the corpus
- `get_all_tissue_ids_in_corpus(snapshot)` - Gets all tissue IDs present in the corpus
- `setup_retry_session()` - Creates a requests session with retry logic

---

### 2. `backend.common.marker_genes`

#### 2.1 `backend.common.marker_genes.computational_markers`
**Used in:**
- `backend/cellguide/pipeline/computational_marker_genes/__init__.py`

**Classes used:**
- `MarkerGenesCalculator` - Calculates computational marker genes from census data
  - `get_computational_marker_genes()` - Computes marker genes based on groupby terms

#### 2.2 `backend.common.marker_genes.constants`
**Used in:**
- `backend/cellguide/pipeline/computational_marker_genes/__init__.py`

**Constants used:**
- `MARKER_SCORE_THRESHOLD` - Threshold for filtering marker genes by score

#### 2.3 `backend.common.marker_genes.marker_gene_files.gene_metadata`
**Used in:**
- `backend/cellguide/pipeline/canonical_marker_genes/canonical_markers.py`

**Functions used:**
- `get_gene_id_to_name_and_symbol()` - Gets mapping from gene IDs to names and symbols

---

### 3. `backend.common.utils`

#### 3.1 `backend.common.utils.dataclass`
**Used in:**
- `backend/cellguide/pipeline/utils.py`

**Functions used:**
- `convert_dataclass_to_dict(data)` - Converts dataclass instances to dictionaries for JSON serialization

#### 3.2 `backend.common.utils.cloudfront`
**Used in:**
- `backend/cellguide/pipeline/__init__.py`

**Functions used:**
- `create_invalidation_for_cellguide_data()` - Invalidates CloudFront cache for CellGuide data after pipeline runs

#### 3.3 `backend.common.utils.result_notification`
**Used in:**
- `backend/cellguide/pipeline/__main__.py`

**Functions used:**
- `format_failed_batch_issue_slack_alert(message)` - Formats failure messages for Slack
- `gen_cg_pipeline_failure_message(message)` - Generates CellGuide pipeline failure message
- `gen_cg_pipeline_success_message(output_path, description_output_path)` - Generates success message
- `notify_slack(message)` - Sends notification to Slack

---

### 4. `backend.common.doi`
**Used in:**
- `backend/cellguide/pipeline/canonical_marker_genes/canonical_markers.py`

**Functions used:**
- `clean_doi(doi)` - Cleans and normalizes DOI strings

---

### 5. `backend.common.providers`

#### 5.1 `backend.common.providers.crossref_provider`
**Used in:**
- `backend/cellguide/pipeline/canonical_marker_genes/canonical_markers.py`

**Classes used:**
- `CrossrefProvider` - Provider for fetching publication metadata from Crossref API
  - `get_title_and_citation_from_doi(doi)` - Fetches title and citation for a DOI

---

## Usage by Pipeline Module

### Main Pipeline (`__init__.py`, `__main__.py`)
- `backend.common.utils.cloudfront.create_invalidation_for_cellguide_data`
- `backend.common.utils.result_notification.*` (4 functions)

### Utils (`utils.py`)
- `backend.common.utils.dataclass.convert_dataclass_to_dict`

### Metadata Pipeline (`metadata/`)
- `backend.common.census_cube.data.snapshot.load_snapshot`
- `backend.common.census_cube.utils.get_all_cell_type_ids_in_corpus`
- `backend.common.census_cube.utils.get_all_tissue_ids_in_corpus`
- `backend.common.census_cube.data.ontology_labels.*` (4 functions)

### Ontology Tree Pipeline (`ontology_tree/`)
- `backend.common.census_cube.data.snapshot.load_snapshot`
- `backend.common.census_cube.data.ontology_labels.ontology_term_label`
- `backend.common.census_cube.utils.*` (7 functions)

### Source Collections Pipeline (`source_collections/`)
- `backend.common.census_cube.data.snapshot.load_snapshot`
- `backend.common.census_cube.data.ontology_labels.ontology_term_label`
- `backend.common.census_cube.utils.descendants`
- `backend.common.census_cube.utils.get_all_cell_type_ids_in_corpus`

### Computational Marker Genes Pipeline (`computational_marker_genes/`)
- `backend.common.census_cube.data.snapshot.load_snapshot`
- `backend.common.marker_genes.computational_markers.MarkerGenesCalculator`
- `backend.common.marker_genes.constants.MARKER_SCORE_THRESHOLD`

### Canonical Marker Genes Pipeline (`canonical_marker_genes/`)
- `backend.common.census_cube.data.snapshot.load_snapshot`
- `backend.common.census_cube.data.ontology_labels.ontology_term_label`
- `backend.common.census_cube.utils.setup_retry_session`
- `backend.common.census_cube.utils.get_all_cell_type_ids_in_corpus`
- `backend.common.doi.clean_doi`
- `backend.common.marker_genes.marker_gene_files.gene_metadata.get_gene_id_to_name_and_symbol`
- `backend.common.providers.crossref_provider.CrossrefProvider`

### GPT Descriptions Pipeline (`gpt_descriptions/`)
- `backend.common.census_cube.data.snapshot.load_snapshot`

---

## Key Dependencies

The most heavily used `backend.common` modules are:

1. **`backend.common.census_cube`** - Core dependency for accessing census data and ontology utilities
2. **`backend.common.marker_genes`** - Used for computational and canonical marker gene processing
3. **`backend.common.utils`** - Utility functions for data serialization, notifications, and cloud operations

The pipeline relies heavily on the census cube snapshot for accessing cell counts, ontology data, and primary filter dimensions, which are used across multiple sub-pipelines.
