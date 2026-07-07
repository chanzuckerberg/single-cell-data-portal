# Shared Dependency Migration Options

## Problem Statement

The `backend/common/` module is shared between:
- **APIs** (DE, WMG, main API server) - likely staying in this repo
- **Pipelines** (CellGuide, WMG pipeline) - potentially moving to separate repo(s)

During migration to a new architecture, these components may not exist in the same repository, requiring a strategy to manage shared dependencies.

---

## Current Dependency Analysis

### What Pipelines Need from Common

| Module | CellGuide Pipeline | WMG Pipeline |
|--------|:------------------:|:------------:|
| `census_cube.data.snapshot` | ✅ | ✅ |
| `census_cube.data.ontology_labels` | ✅ | ✅ |
| `census_cube.utils` | ✅ | ✅ |
| `marker_genes.computational_markers` | ✅ | ✅ |
| `marker_genes.blacklist` | - | ✅ |
| `utils.result_notification` | ✅ | ✅ |
| `utils.cloudfront` | ✅ | - |
| `utils.secret_config` | ✅ | - |
| `providers.crossref_provider` | ✅ | - |

### What APIs Need (for comparison)

| Module | DE API | WMG API |
|--------|:------:|:-------:|
| `server.config` | ✅ | ✅ |
| `census_cube.data.*` | ✅ | ✅ |
| `marker_genes.blacklist` | ✅ | ✅ |

---

## Option 1: Published Python Package

**Approach:** Extract `backend/common/` into a standalone Python package published to a private PyPI registry (or public if appropriate).

```
# New package structure
cellxgene-common/
├── src/
│   └── cellxgene_common/
│       ├── census_cube/
│       ├── marker_genes/
│       ├── server/
│       └── utils/
├── pyproject.toml
└── README.md
```

**Usage in consuming repos:**
```python
# pyproject.toml or requirements.txt
cellxgene-common>=1.0.0

# Code
from cellxgene_common.census_cube.data import snapshot
from cellxgene_common.marker_genes import MarkerGenesCalculator
```

### Pros
- Clean versioning with semantic versioning
- Standard Python packaging workflow
- Independent release cycle
- Clear API boundaries enforced by package structure
- Works with any CI/CD system
- Dependency pinning prevents unexpected breaks

### Cons
- Requires private PyPI infrastructure (AWS CodeArtifact, Artifactory, etc.)
- Version coordination overhead between repos
- Local development requires editable installs or frequent publishes
- Breaking changes require coordinated releases

### Implementation Steps
1. Set up private PyPI (AWS CodeArtifact recommended for AWS-native stack)
2. Extract common code into new package repo
3. Add CI/CD to publish on merge to main
4. Update consuming repos to install from private registry
5. Establish versioning policy (semver recommended)

### Best For
- Long-term maintainability
- Multiple teams working on different repos
- Need for stable, versioned dependencies

---

## Option 2: Git Submodule

**Approach:** Keep common code in its own repo, reference it as a git submodule in consuming repos.

```
# In pipeline repo
pipeline-repo/
├── .gitmodules
├── shared/
│   └── common/  # git submodule → common-lib repo
├── pipeline/
│   └── cellguide/
└── pyproject.toml
```

```ini
# .gitmodules
[submodule "shared/common"]
    path = shared/common
    url = git@github.com:org/cellxgene-common.git
    branch = main
```

**Usage:**
```python
# Add to Python path or use relative imports
import sys
sys.path.insert(0, 'shared/common')

from census_cube.data import snapshot
```

### Pros
- No package registry infrastructure needed
- Always working with source code (easy debugging)
- Can pin to specific commits
- Familiar git workflow

### Cons
- Submodules are notoriously confusing for developers
- Easy to get into inconsistent states
- Requires manual updates (`git submodule update`)
- CI/CD needs special handling for submodule checkout
- No true versioning, just commit pinning

### Implementation Steps
1. Extract common code to separate repo
2. Add as submodule in consuming repos
3. Configure CI/CD to checkout submodules
4. Add path manipulation for Python imports

### Best For
- Small teams comfortable with git submodules
- Rapid iteration where package publishing is too slow
- Temporary solution during migration

---

## Option 3: Monorepo with Published Artifacts

**Approach:** Keep everything in the monorepo but publish common code as a package from within the monorepo.

```
single-cell-data-portal/
├── backend/
│   ├── common/           # Published as package
│   ├── api_server/       # Imports common locally OR as package
│   ├── de/
│   ├── wmg/
│   └── cellguide/
├── pipelines/            # Could be separate CI/CD
│   └── ...
└── packages/
    └── cellxgene-common/
        └── pyproject.toml  # Points to backend/common
```

### Pros
- Single source of truth
- Local development is simple (relative imports)
- Atomic commits across common + consumers
- Can publish for external consumers

### Cons
- Doesn't solve the "separate repo" requirement
- Complex CI/CD (conditional publishing)
- Version management within monorepo is tricky

### Best For
- If pipelines can stay in monorepo with separate deployment
- Hybrid approach during migration

---

## Option 4: Multiple Domain-Specific Packages

**Approach:** Split `common/` into multiple focused packages by domain, allowing repos to depend only on what they need.

```
Packages:
├── cellxgene-census-cube    # census_cube module
├── cellxgene-marker-genes   # marker_genes module
├── cellxgene-server-utils   # server + API utilities
└── cellxgene-pipeline-utils # pipeline-specific utils
```

**Dependency graph:**
```
                    ┌─────────────────────┐
                    │ cellxgene-census-cube│
                    └──────────┬──────────┘
                               │
            ┌──────────────────┼──────────────────┐
            │                  │                  │
            ▼                  ▼                  ▼
┌───────────────────┐ ┌───────────────┐ ┌─────────────────┐
│cellxgene-marker-  │ │cellxgene-     │ │cellxgene-       │
│genes              │ │server-utils   │ │pipeline-utils   │
└─────────┬─────────┘ └───────┬───────┘ └────────┬────────┘
          │                   │                  │
          │                   ▼                  │
          │           ┌──────────────┐           │
          │           │   DE API     │           │
          │           │   WMG API    │           │
          │           └──────────────┘           │
          │                                      │
          └──────────────┬───────────────────────┘
                         ▼
               ┌──────────────────┐
               │ CellGuide Pipeline│
               │ WMG Pipeline      │
               └──────────────────┘
```

### Pros
- Fine-grained dependencies (install only what you need)
- Clearer ownership boundaries
- Independent versioning per domain
- Smaller package sizes

### Cons
- More packages to maintain
- Potential for circular dependency issues
- More complex release coordination
- Inter-package version compatibility matrix

### Suggested Package Split

| Package | Contents | Used By |
|---------|----------|---------|
| `cellxgene-census-cube` | `census_cube/` | All |
| `cellxgene-marker-genes` | `marker_genes/` | CellGuide, WMG |
| `cellxgene-server` | `server/`, `corpora_config` | APIs only |
| `cellxgene-utils` | `utils/`, `doi`, `providers` | All |

### Best For
- Large organizations with clear domain ownership
- When different components have different release cadences
- Minimizing dependency footprint

---

## Option 5: Service Abstraction (API-ify Shared Logic)

**Approach:** Convert shared functionality into internal microservices rather than shared libraries.

```
Services:
├── census-cube-service     # Exposes census data via API
├── marker-genes-service    # Computes marker genes on demand
└── ontology-service        # Ontology lookups
```

**Example:**
```python
# Instead of:
from cellxgene_common.census_cube.data import snapshot
data = snapshot.load_snapshot()

# Use:
import httpx
response = httpx.get("http://census-cube-service/snapshot")
data = response.json()
```

### Pros
- No shared code dependencies at all
- Services can be scaled independently
- Language-agnostic (pipelines could be in different languages)
- Clear API contracts

### Cons
- Significant architectural change
- Network latency for what was local calls
- More infrastructure to maintain
- Overkill for simple utility functions
- Not suitable for heavy computation (marker genes)

### What Could Be Service-ified

| Functionality | Good Candidate? | Reason |
|---------------|-----------------|--------|
| Ontology lookups | ✅ Yes | Simple, cacheable |
| Snapshot loading | ⚠️ Maybe | Large data transfer |
| Marker gene calculation | ❌ No | Computation-heavy, needs local data |
| Slack notifications | ✅ Yes | Simple, async-friendly |

### Best For
- When shared logic is already stateless and RPC-like
- Multi-language environments
- When you want to enforce API boundaries strictly

---

## Option 6: Copy with Automated Sync

**Approach:** Copy shared code into consuming repos with automated synchronization.

```yaml
# .github/workflows/sync-common.yml
name: Sync Common Code
on:
  repository_dispatch:
    types: [common-updated]
jobs:
  sync:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - name: Copy from common repo
        run: |
          git clone https://github.com/org/common-lib.git /tmp/common
          cp -r /tmp/common/src/* ./shared/
      - name: Create PR
        uses: peter-evans/create-pull-request@v5
```

### Pros
- No submodule complexity
- No package registry needed
- Full code visibility in each repo
- Can diverge if needed (escape hatch)

### Cons
- Code duplication
- Sync can fail or drift
- Merge conflicts when syncing
- No clear versioning

### Best For
- Temporary migration state
- When you expect repos to eventually diverge
- Small amount of shared code

---

## Recommendation

### For This Codebase: Option 1 (Published Package) + Option 4 (Domain Split)

**Rationale:**
1. The shared code is substantial (59 files) and well-structured by domain
2. Clear separation already exists (`census_cube`, `marker_genes`, `server`, `utils`)
3. APIs and pipelines have different dependency needs
4. AWS-native stack makes CodeArtifact a natural fit

### Suggested Implementation

**Phase 1: Single Package**
```
cellxgene-backend-common
├── census_cube/
├── marker_genes/
├── server/
├── utils/
└── providers/
```

**Phase 2: Split if Needed**
```
cellxgene-census-cube     # Core data (used by all)
cellxgene-marker-genes    # Depends on census-cube
cellxgene-backend-utils   # Standalone utilities
```

### Migration Path

```
Week 1-2: Setup
├── Create private PyPI (AWS CodeArtifact)
├── Extract common to new repo
├── Set up CI/CD for publishing
└── Version: 0.1.0

Week 3-4: Integration
├── Update main repo to use package
├── Validate all tests pass
├── Document installation for developers
└── Version: 0.2.0 (if changes needed)

Week 5+: Pipeline Migration
├── New pipeline repo installs package
├── Iterate on package API as needed
└── Version: 1.0.0 when stable
```

### Local Development Setup

```bash
# For package development
cd cellxgene-common
pip install -e ".[dev]"

# For consuming repo (editable install of local common)
cd single-cell-data-portal
pip install -e ../cellxgene-common
```

---

## Comparison Matrix

| Criteria | Package | Submodule | Monorepo | Multi-Pkg | Service | Copy/Sync |
|----------|:-------:|:---------:|:--------:|:---------:|:-------:|:---------:|
| Setup complexity | Medium | Low | Low | High | Very High | Low |
| Maintenance burden | Low | Medium | Low | High | High | High |
| Version control | ✅ | ⚠️ | ❌ | ✅ | ✅ | ❌ |
| Local dev experience | Good | Poor | Great | Good | Poor | Good |
| CI/CD complexity | Medium | Medium | Low | High | Very High | Medium |
| Scales to many repos | ✅ | ⚠️ | ❌ | ✅ | ✅ | ❌ |
| Breaking change safety | ✅ | ❌ | ❌ | ✅ | ✅ | ❌ |

**Legend:** ✅ = Good, ⚠️ = Caution, ❌ = Poor

---

## Quick Start: AWS CodeArtifact Setup

```bash
# Create repository
aws codeartifact create-repository \
  --domain cellxgene \
  --repository python-packages

# Configure pip
aws codeartifact login --tool pip \
  --domain cellxgene \
  --repository python-packages

# Publish package
cd cellxgene-common
pip install build twine
python -m build
twine upload --repository codeartifact dist/*
```

```toml
# pyproject.toml for consuming repo
[project]
dependencies = [
    "cellxgene-backend-common>=1.0.0",
]

[tool.pip]
extra-index-url = "https://cellxgene-ACCOUNT.d.codeartifact.REGION.amazonaws.com/pypi/python-packages/simple/"
```
