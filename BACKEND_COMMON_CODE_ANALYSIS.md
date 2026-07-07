# Backend Common Code Analysis

## Overview

This document maps the shared common code in `backend/common/` and how it's used across the CellGuide, Differential Expression (DE), and Where's My Gene (WMG) modules.

**Directories Analyzed:**
- `backend/common/` - Shared utilities and infrastructure (59 Python files)
- `backend/cellguide/` - CellGuide pipeline and services (43 files)
- `backend/de/` - Differential Expression API (3 files)
- `backend/wmg/` - Where's My Gene API and pipeline (29 files)

---

## 1. Common Module Structure

```
backend/common/
├── census_cube/          # TileDB cube data management (18 files)
│   ├── data/
│   │   ├── snapshot.py       # CensusCubeSnapshot dataclass
│   │   ├── query.py          # CensusCubeQuery class
│   │   ├── criteria.py       # Query criteria models
│   │   ├── ontology_labels.py # Ontology term lookups
│   │   ├── constants.py      # Assays, schema version
│   │   ├── tiledb.py         # TileDB context
│   │   └── schemas/          # Cube schema definitions
│   ├── utils.py              # Filter/dimension utilities
│   └── config.py             # CensusCubeConfig
├── marker_genes/         # Marker gene calculations (8 files)
│   ├── computational_markers.py  # MarkerGenesCalculator
│   ├── types.py              # ComputationalMarkerGenes dataclass
│   ├── constants.py          # MARKER_SCORE_THRESHOLD
│   ├── utils.py              # Calculation utilities
│   └── marker_gene_files/    # Blacklists, gene metadata
├── server/               # Flask/API configuration (4 files)
│   ├── config.py             # create_api_app factory
│   ├── logger.py             # Logging configuration
│   ├── request_id.py         # Request ID management
│   └── datadog.py            # Tracing initialization
├── utils/                # General utilities (26 files)
│   ├── secret_config.py      # AWS Secrets Manager base
│   ├── http_exceptions.py    # Custom HTTP exceptions
│   ├── result_notification.py # Slack notifications
│   ├── cloudfront.py         # CDN invalidation
│   ├── dataclass.py          # Dataclass utilities
│   ├── math_utils.py         # Math constants
│   └── ...                   # Many more utilities
├── providers/            # Third-party integrations
│   └── crossref_provider.py  # CrossRef API
├── corpora_config.py     # Main app configuration
├── auth0_manager.py      # Auth0 JWT validation
├── authorizer.py         # Authorization utilities
├── citation.py           # Citation formatting
├── doi.py                # DOI handling
├── feature_flag.py       # Feature flags
├── constants.py          # Global constants
└── logging_config.py     # Logging constants
```

---

## 2. Core Modules Detail

### A. Census Cube (`census_cube/`)

The core data infrastructure for expression summaries, cell counts, and marker genes.

| File | Purpose | Key Exports |
|------|---------|-------------|
| `data/snapshot.py` | Census data artifacts | `CensusCubeSnapshot`, `load_snapshot()` |
| `data/query.py` | Query execution | `CensusCubeQuery`, `CensusCubeQueryParams` |
| `data/criteria.py` | Query filtering | `BaseQueryCriteria`, `CensusCubeQueryCriteria` |
| `data/ontology_labels.py` | Ontology lookups | `ontology_term_label()`, `ontology_term_description()` |
| `data/constants.py` | Configuration | `INCLUDED_ASSAYS`, `CENSUS_CUBE_DATA_SCHEMA_VERSION` |
| `utils.py` | Utilities | `ancestors()`, `descendants()`, `find_all_dim_option_values()` |

### B. Marker Genes (`marker_genes/`)

Computational marker gene calculations for cell types.

| File | Purpose | Key Exports |
|------|---------|-------------|
| `computational_markers.py` | Main calculator | `MarkerGenesCalculator` |
| `types.py` | Data types | `ComputationalMarkerGenes` |
| `constants.py` | Thresholds | `MARKER_SCORE_THRESHOLD`, `PIPELINE_NUM_CPUS` |
| `utils.py` | Calculations | `calculate_specificity_excluding_nans()`, `calculate_cohens_d()` |
| `marker_gene_files/blacklist.py` | Exclusions | Ribosomal, mitochondrial gene lists |

### C. Server (`server/`)

Flask/Connexion server configuration shared by API services.

| File | Purpose | Key Exports |
|------|---------|-------------|
| `config.py` | App factory | `create_api_app()`, `configure_flask_app()` |
| `logger.py` | Logging | `configure_logging()` |
| `request_id.py` | Request tracking | Request ID generation |
| `datadog.py` | Tracing | Datadog initialization |

### D. Utils (`utils/`)

General utilities shared across all backend services.

| File | Purpose | Key Exports |
|------|---------|-------------|
| `secret_config.py` | Config base | `SecretConfig` (AWS Secrets Manager) |
| `http_exceptions.py` | HTTP errors | `UnauthorizedError`, `ForbiddenHTTPException`, etc. |
| `result_notification.py` | Notifications | `notify_slack()`, `upload_to_slack()` |
| `cloudfront.py` | CDN | CloudFront invalidation |
| `dataclass.py` | Conversion | `convert_dataclass_to_dict()` |
| `math_utils.py` | Constants | `GB`, `MB` |

---

## 3. Dependency Matrix

### By Common Module

| Common Module | CellGuide | DE | WMG | Description |
|---------------|:---------:|:--:|:---:|-------------|
| `census_cube.data.snapshot` | ✅ | ✅ | ✅ | Core data loading |
| `census_cube.data.query` | - | ✅ | ✅ | Query execution |
| `census_cube.data.criteria` | - | ✅ | ✅ | Query criteria |
| `census_cube.data.ontology_labels` | ✅ | ✅ | ✅ | Term lookups |
| `census_cube.data.tiledb` | - | - | ✅ | TileDB context |
| `census_cube.data.schemas` | - | ✅ | ✅ | Cube schemas |
| `census_cube.utils` | ✅ | ✅ | ✅ | Filter/dimension utilities |
| `marker_genes.computational_markers` | ✅ | - | ✅ | Marker calculations |
| `marker_genes.constants` | ✅ | - | - | Score threshold |
| `marker_genes.blacklist` | - | ✅ | ✅ | Gene exclusion |
| `server.config` | - | ✅ | ✅ | App factory |
| `utils.cloudfront` | ✅ | - | - | CDN invalidation |
| `utils.result_notification` | ✅ | - | ✅ | Slack notifications |
| `utils.secret_config` | ✅ | - | - | Configuration |
| `utils.http_exceptions` | ✅ | - | - | HTTP errors |
| `utils.dataclass` | ✅ | - | - | Data conversion |
| `utils.math_utils` | - | - | ✅ | Constants |
| `doi` | ✅ | - | - | DOI cleaning |
| `providers.crossref_provider` | ✅ | - | - | Citation lookup |

### By Application

#### CellGuide Imports
```python
from backend.common.census_cube.data import snapshot, ontology_labels
from backend.common.census_cube.utils import (
    get_all_cell_type_ids_in_corpus,
    get_all_tissue_ids_in_corpus,
    descendants,
    setup_retry_session
)
from backend.common.marker_genes import MarkerGenesCalculator, MARKER_SCORE_THRESHOLD
from backend.common.marker_genes.marker_gene_files import gene_metadata
from backend.common.utils import cloudfront, dataclass, http_exceptions, result_notification, secret_config
from backend.common.doi import clean_doi
from backend.common.providers import CrossrefProvider
```

#### DE Imports
```python
from backend.common.server.config import create_api_app
from backend.common.census_cube.data import criteria, ontology_labels, query, snapshot, schemas
from backend.common.census_cube.utils import ancestors, descendants
from backend.common.marker_genes import marker_gene_blacklist
```

#### WMG Imports
```python
from backend.common.server.config import create_api_app
from backend.common.census_cube.data import constants, criteria, ontology_labels, query, snapshot, schemas, tiledb
from backend.common.census_cube.utils import (
    ancestors,
    descendants,
    are_cell_types_not_redundant_nodes,
    build_filter_relationships
)
from backend.common.marker_genes import MarkerGenesCalculator, computational_markers
from backend.common.utils import result_notification, math_utils
```

---

## 4. Usage Patterns

### Universally Shared (All 3 Apps)

| Module | Function |
|--------|----------|
| `census_cube.data.snapshot` | Loading census cube data artifacts |
| `census_cube.data.ontology_labels` | Looking up ontology term labels/descriptions |
| `census_cube.utils` | Cell type hierarchy traversal (`ancestors`, `descendants`) |

### API Services Only (DE + WMG)

| Module | Function |
|--------|----------|
| `server.config` | Flask app factory (`create_api_app()`) |
| `census_cube.data.query` | Query execution against TileDB |
| `census_cube.data.criteria` | Query parameter models |
| `census_cube.data.schemas` | Cube schema definitions |
| `marker_genes.blacklist` | Gene exclusion lists |

### Pipeline Only (CellGuide + WMG)

| Module | Function |
|--------|----------|
| `marker_genes.computational_markers` | Marker gene calculations |
| `utils.result_notification` | Slack notifications for pipeline status |

### CellGuide Exclusive

| Module | Function |
|--------|----------|
| `utils.cloudfront` | CDN invalidation after data updates |
| `providers.crossref_provider` | Citation lookups from CrossRef |
| `doi` | DOI string cleaning |
| `utils.secret_config` | Direct configuration access |

---

## 5. Dependency Hierarchy

### Tier 1: Core Infrastructure (Used by All)
```
census_cube.data.snapshot
census_cube.data.ontology_labels
census_cube.utils
```

### Tier 2: API Infrastructure (Used by DE + WMG)
```
server.config
census_cube.data.query
census_cube.data.criteria
census_cube.data.schemas
```

### Tier 3: Computation (Used by CellGuide + WMG)
```
marker_genes.computational_markers
marker_genes.types
marker_genes.utils
```

### Tier 4: Pipeline Support
```
utils.result_notification (CellGuide + WMG)
utils.cloudfront (CellGuide only)
providers.crossref_provider (CellGuide only)
```

---

## 6. Layers Integration

The `backend/layers/` directory provides additional shared functionality:

| App | Layers Usage |
|-----|--------------|
| CellGuide | Yes - `backend.layers` for upload actions |
| DE | No |
| WMG | No |

Layers provides: auth, business logic, persistence, and processing utilities primarily used by the main API server and CellGuide.

---

## 7. Summary Statistics

| Metric | Count |
|--------|-------|
| Total common Python files | 59 |
| CellGuide Python files | 43 |
| DE Python files | 3 |
| WMG Python files | 29 |

### Common Modules by Usage Count

| Modules | Used By | Count |
|---------|---------|-------|
| `census_cube.data.snapshot` | All | 3 |
| `census_cube.data.ontology_labels` | All | 3 |
| `census_cube.utils` | All | 3 |
| `server.config` | DE, WMG | 2 |
| `census_cube.data.query` | DE, WMG | 2 |
| `census_cube.data.criteria` | DE, WMG | 2 |
| `marker_genes.computational_markers` | CellGuide, WMG | 2 |
| `marker_genes.blacklist` | DE, WMG | 2 |
| `utils.result_notification` | CellGuide, WMG | 2 |

---

## 8. Architecture Diagram

```
                    ┌─────────────────────────────────────────┐
                    │            backend/common/              │
                    └─────────────────────────────────────────┘
                                       │
         ┌─────────────────────────────┼─────────────────────────────┐
         │                             │                             │
         ▼                             ▼                             ▼
┌─────────────────┐         ┌─────────────────┐         ┌─────────────────┐
│   census_cube   │         │  marker_genes   │         │     server      │
│                 │         │                 │         │                 │
│ • snapshot      │         │ • calculator    │         │ • create_api_app│
│ • query         │         │ • blacklist     │         │ • logging       │
│ • criteria      │         │ • types         │         │ • request_id    │
│ • ontology      │         │ • constants     │         │                 │
│ • schemas       │         │                 │         │                 │
└────────┬────────┘         └────────┬────────┘         └────────┬────────┘
         │                           │                           │
         └───────────┬───────────────┴───────────────────────────┘
                     │
    ┌────────────────┼────────────────┬────────────────┐
    │                │                │                │
    ▼                ▼                ▼                ▼
┌────────┐    ┌────────────┐    ┌──────────┐    ┌──────────┐
│  utils │    │ providers  │    │   doi    │    │ corpora  │
│        │    │            │    │          │    │ _config  │
│• secret│    │• crossref  │    │• clean   │    │          │
│• http  │    │            │    │• normalize│   │• secrets │
│• slack │    │            │    │          │    │• db      │
│• cloud │    │            │    │          │    │          │
└────────┘    └────────────┘    └──────────┘    └──────────┘
                     │
    ┌────────────────┼────────────────┬────────────────┐
    │                │                │                │
    ▼                ▼                ▼                ▼
┌──────────┐   ┌──────────┐    ┌──────────┐    ┌──────────┐
│CellGuide │   │    DE    │    │   WMG    │    │API Server│
│(pipeline)│   │  (API)   │    │(API+pipe)│    │  (main)  │
└──────────┘   └──────────┘    └──────────┘    └──────────┘
```

---

## 9. Key Observations

1. **Census Cube is Central**: The `census_cube` module is the most critical shared infrastructure, used by all three applications for data access.

2. **Clear API vs Pipeline Split**: DE and WMG use `server.config` for API setup, while CellGuide focuses on pipeline utilities.

3. **Marker Genes Shared for Computation**: Both CellGuide and WMG perform marker gene calculations using the shared `MarkerGenesCalculator`.

4. **CellGuide is Most Diverse**: CellGuide uses the widest variety of common modules, including exclusive use of CloudFront invalidation and CrossRef integration.

5. **DE is Most Focused**: With only 3 Python files, DE has the narrowest dependency footprint, focusing purely on query execution.

6. **Ontology Labels Universal**: All three apps need to resolve ontology term IDs to human-readable labels.
