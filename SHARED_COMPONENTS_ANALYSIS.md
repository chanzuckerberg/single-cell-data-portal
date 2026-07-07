# Shared Components Analysis: DifferentialExpression & WheresMyGeneV2

## Overview

This document analyzes shared components between the DifferentialExpression (DE) and WheresMyGeneV2 (WMG) applications, identifying truly shared components, duplicated patterns, and opportunities for consolidation.

**Entry Points Analyzed:**
- `frontend/src/views/DifferentialExpression/components/Main/index.tsx`
- `frontend/src/views/WheresMyGeneV2/components/Main/index.tsx`

---

## 1. Directly Shared Components from `src/components/`

### BottomBanner
**Location:** `src/components/BottomBanner`

| App | Usage |
|-----|-------|
| DifferentialExpression | Main component |
| WheresMyGeneV2 | Main component |
| CellGuide | Also uses this |

**Purpose:** Newsletter signup & feedback survey banner displayed at bottom of page.

### GeneInfoSideBar
**Location:** `src/components/GeneInfoSideBar`

| App | Usage |
|-----|-------|
| WheresMyGeneV2 | Main component (right sidebar) |
| CellGuide | CellGuideCard component |

**Purpose:** Displays gene information fetched from NCBI.

### Notification
**Location:** `src/components/Notification`

| App | Usage |
|-----|-------|
| WheresMyGeneV2 | Main, CitationButton, ShareButton |

**Purpose:** Toast/notification system using WMG's StateContext.

### SideBar
**Location:** `src/components/common/SideBar`

| App | Usage |
|-----|-------|
| WheresMyGeneV2 | Filters panel |
| Datasets | Filter sidebar |
| Collections | Filter sidebar |

**Purpose:** Collapsible sidebar with toggle button.

---

## 2. Shared Constants & Utilities

### Analytics
**Location:** `src/common/analytics`

Both applications use:
- `track` function from `src/common/analytics`
- `EVENTS` from `src/common/analytics/events`

### Constants
| Constant | Location | Used By |
|----------|----------|---------|
| `BANNER_FEEDBACK_SURVEY_LINK` | `src/common/constants/airtableLinks` | DE, WMG |
| `EMPTY_ARRAY`, `noop` | `src/common/constants/utils` | Both |
| `ROUTES` | `src/common/constants/routes` | Both |

### Theme Variables
**Location:** `src/common/theme`

Both use color and spacing variables:
- `primary400`, `gray500`, `gray400`
- Various spacing constants

### Layout Constants
**Location:** `src/components/Layout/style`

- `CONTENT_WRAPPER_LEFT_RIGHT_PADDING_PX`
- `HEADER_HEIGHT_PX`
- `FOOTER_HEIGHT_PX`
- `SidebarMainWrapper`

---

## 3. Duplicated Components (Candidates for Consolidation)

### Loader Components

| App | Location |
|-----|----------|
| DifferentialExpression | `components/Main/components/Loader` |
| WheresMyGeneV2 | `components/Loader` |

Both wrap `LoadingIndicator` from `@czi-sds/components` with similar implementations.

**Recommendation:** Create unified `src/components/Loader`.

### Organism Type Definition

Both apps define identical interfaces:

```typescript
interface Organism {
  id: string;
  name: string;
}
```

| App | Location |
|-----|----------|
| DifferentialExpression | `common/types.ts` |
| WheresMyGeneV2 | `common/types.ts` |

**Recommendation:** Move to `src/common/types/organism.ts`.

### Organism Dropdown Components

| App | Location | Implementation |
|-----|----------|----------------|
| DifferentialExpression | `components/Main/components/Organism` | Simple dropdown |
| WheresMyGeneV2 | `components/Filters/components/Organism` | Dropdown + Dialog |

Both follow similar logic but with different UI patterns.

**Recommendation:** Create `src/components/OrganismSelector` with variants.

---

## 4. Similar Patterns with Different Implementations

### Filter Systems

| App | Approach |
|-----|----------|
| DifferentialExpression | Custom `FilterDropdown` with `FilterGroup` logic |
| WheresMyGeneV2 | `StyledComplexFilter` from CZI-SDS with Organism, Compare, Sort, ColorScale subcomponents |

### State Management

Both use Redux-like patterns:

| App | Store Location |
|-----|----------------|
| DifferentialExpression | `common/store/` (reducer, actions) |
| WheresMyGeneV2 | `common/store/` (reducer, actions, selectors, StateContext) |

**Shared Pattern:** Both use `useConnect` hooks to separate business logic from components.

---

## 5. Cross-App Dependencies

### Important Finding

DifferentialExpression imports from WheresMyGeneV2:
- `StyledSidebarDrawer` from `WheresMyGeneV2/components/Main/style` (used in `DeResults`)

**Recommendation:** Move `StyledSidebarDrawer` to `src/components/common/`.

---

## 6. Components Used by Other Views

| Component | DE | WMG | CellGuide | Datasets | Collections |
|-----------|:--:|:---:|:---------:|:--------:|:-----------:|
| BottomBanner | X | X | X | - | - |
| GeneInfoSideBar | - | X | X | - | - |
| SideBar | - | X | - | X | X |
| Notification | - | X | - | - | - |
| Table (from CellGuide) | X | - | X | - | - |
| HelpTooltip (from CellGuide) | X | - | X | - | - |

---

## 7. Summary Table

| Component/Feature | DE | WMG | Status |
|-------------------|:--:|:---:|--------|
| BottomBanner | X | X | **Actively shared** |
| GeneInfoSideBar | - | X | Partial (via CellGuide) |
| Notification | - | X | WMG-specific |
| SideBar | - | X | Different approaches |
| Loader | X | X | **Duplicated - consolidate** |
| Organism dropdown | X | X | **Similar - consolidate** |
| Organism type | X | X | **Duplicated - consolidate** |
| Filters system | X | X | Different implementations |
| Table component | X | - | Uses CellGuide's Table |
| Analytics | X | X | **Actively shared** |
| Theme/Colors | X | X | **Actively shared** |
| Routes | X | X | **Actively shared** |

---

## 8. Recommendations

### High Priority (Quick Wins)

1. **Move Organism type to shared location**
   - Create `src/common/types/organism.ts`
   - Import in both DE and WMG

2. **Consolidate Loader components**
   - Create `src/components/Loader`
   - Both currently wrap `LoadingIndicator` from SDS

3. **Move StyledSidebarDrawer to shared**
   - Currently DE imports from WMG's style file
   - Move to `src/components/common/SidebarDrawer`

### Medium Priority

4. **Standardize dropdown styling**
   - Create unified styled dropdown variants in `src/components/common/Dropdown`

5. **Extract organism selection logic**
   - Create `src/components/OrganismSelector` with basic and advanced variants

### Lower Priority

6. **Unify notification system**
   - DE could benefit from WMG's notification system
   - Requires StateContext provider wrapper

7. **Consider shared filter architecture**
   - Long-term: Create abstract filter component
   - Both need similar filter-to-results workflows

---

## 9. Conclusion

| Category | Count |
|----------|-------|
| Truly Shared Components | 3 (BottomBanner, analytics, theme constants) |
| Duplicated/Similar | 3 (Loader, Organism types, Organism component) |
| Different Implementations | 2 (Filters systems, state management) |

The applications share more **architectural patterns** and **utilities** than actual components. The main opportunities for consolidation are:

1. Type definitions (Organism)
2. Loading states (Loader component)
3. Drawer/sidebar UI components
4. Filter building blocks

Both applications are well-structured but could benefit from more shared infrastructure in the middle layer between app-specific code and the base component library.
