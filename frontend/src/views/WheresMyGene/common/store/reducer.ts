import isEqual from "lodash/isEqual";
import { CompareId, X_AXIS_CHART_HEIGHT_PX } from "../constants";
import { CellType, SORT_BY } from "../types";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { GENE_SEARCH_BAR_HEIGHT_PX } from "src/views/WheresMyGeneV2/common/constants";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";

export interface PayloadAction<Payload> {
  type: keyof typeof REDUCERS;
  payload: Payload;
}

export interface State {
  genesToDelete: string[];
  selectedGenes: string[];
  selectedOrganismId: string | null;
  selectedTissues?: string[];
  selectedFilters: {
    datasets: string[];
    developmentStages: string[];
    diseases: string[];
    ethnicities: string[];
    publications: string[];
    sexes: string[];
    tissues: string[];
  };

  /**
   * (thuang): BE API response always returns a snapshot ID. When the ID changes,
   * FE needs refresh the queries
   */
  snapshotId: string | null;
  sortBy: { cellTypes: SORT_BY; genes: SORT_BY; scaled: SORT_BY };
  cellInfoCellType: {
    cellType: CellType;
    tissueID: string;
    organismID: string;
  } | null;
  geneInfoGene: string | null;
  compare?: CompareId;
  xAxisHeight: number;
  filteredCellTypes: string[];
  filteredCellTypeIds: string[];
  expandedTissueIds: string[];
}

const EMPTY_FILTERS: State["selectedFilters"] = {
  datasets: [],
  developmentStages: [],
  diseases: [],
  ethnicities: [],
  publications: [],
  sexes: [],
  tissues: [],
};

// (thuang): If you have derived states based on the state, use `useMemo`
// to cache the derived states instead of putting them in the state.
export const INITIAL_STATE: State = {
  cellInfoCellType: null,
  geneInfoGene: null,
  genesToDelete: [],
  selectedFilters: EMPTY_FILTERS,
  selectedGenes: [],
  selectedOrganismId: null,
  selectedTissues: [],
  snapshotId: null,
  sortBy: {
    cellTypes: SORT_BY.CELL_ONTOLOGY,
    genes: SORT_BY.USER_ENTERED,
    scaled: SORT_BY.COLOR_SCALED,
  },
  xAxisHeight: X_AXIS_CHART_HEIGHT_PX + GENE_SEARCH_BAR_HEIGHT_PX,
  filteredCellTypes: [],
  filteredCellTypeIds: [],
  expandedTissueIds: [],
};

export const REDUCERS = {
  addCellInfoCellType,
  addGeneInfoGene,
  clearCellInfoCellType,
  clearGeneInfoGene,
  closeRightSidebar,
  addSelectedGenes,
  deleteSelectedGenes,
  deleteSingleGene,
  deleteAllGenes,
  loadStateFromURL,
  selectCompare,
  resetGenesToDelete,
  selectFilters,
  selectGenes,
  selectGeneInfoFromXAxis,
  selectOrganism,
  selectSortBy,
  selectTissues,
  setSnapshotId,
  toggleGeneToDelete,
  setXAxisHeight,
  setFilteredCellTypes,
  toggleExpandedTissueId,
  autoExpandTissues,
};

export function reducer(state: State, action: PayloadAction<unknown>): State {
  const { type } = action;

  const handler = REDUCERS[type];

  if (!handler) {
    throw new Error(`Unknown action type: ${type}`);
  }

  /**
   * (thuang): Figuring out the typing here is more work than its rewards
   * Ideally we'll use Redux Toolkit's `createSlice` here, but I think that's
   * too heavy for now
   */
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  return handler(state, action as any);
}

/**
 * Individual action type reducers below. Please add their corresponding action
 * creators in `./actions.ts`
 */

function deleteSingleGene(
  state: State,
  action: PayloadAction<string | null>
): State {
  if (!action.payload) {
    return state;
  }

  const { selectedGenes } = state;

  const newSelectedGenes = deleteByItems<State["selectedGenes"][number]>(
    selectedGenes,
    [action.payload]
  );

  return {
    ...state,
    selectedGenes: newSelectedGenes,
    xAxisHeight: INITIAL_STATE.xAxisHeight,
  };
}

function deleteAllGenes(state: State, _: PayloadAction<null>): State {
  return {
    ...state,
    selectedGenes: [],
    xAxisHeight: INITIAL_STATE.xAxisHeight,
  };
}

function deleteSelectedGenes(state: State, _: PayloadAction<null>): State {
  const { genesToDelete } = state;

  if (!genesToDelete.length) {
    return state;
  }

  const { selectedGenes } = state;

  const newSelectedGenes = genesToDelete.length
    ? deleteByItems<State["selectedGenes"][number]>(
        selectedGenes,
        genesToDelete
      )
    : selectedGenes;

  return {
    ...state,
    genesToDelete: [],
    selectedGenes: newSelectedGenes,
  };
}

function selectOrganism(
  state: State,
  action: PayloadAction<string | null>
): State {
  if (state.selectedOrganismId === action.payload) {
    return state;
  }

  return {
    ...state,
    selectedGenes: [],
    selectedOrganismId: action.payload,
    selectedTissues: EMPTY_ARRAY,
    selectedFilters: EMPTY_FILTERS,
    cellInfoCellType: null,
  };
}

function selectGenes(
  state: State,
  action: PayloadAction<State["selectedGenes"]>
): State {
  return {
    ...state,
    genesToDelete: [],
    selectedGenes: action.payload,
  };
}

function selectTissues(
  state: State,
  action: PayloadAction<State["selectedTissues"]>
): State {
  return {
    ...state,
    selectedTissues: action.payload,
  };
}

function selectSortBy(
  state: State,
  action: PayloadAction<Partial<State["sortBy"]>>
): State {
  return {
    ...state,
    sortBy: { ...state.sortBy, ...action.payload },
  };
}

function toggleGeneToDelete(
  state: State,
  action: PayloadAction<string>
): State {
  if (state.genesToDelete.includes(action.payload)) {
    return {
      ...state,
      genesToDelete: deleteByItems<string>(state.genesToDelete, [
        action.payload,
      ]),
    };
  }

  return {
    ...state,
    genesToDelete: [...state.genesToDelete, action.payload],
  };
}

function deleteByItems<Item>(collection: Item[], collectionToDelete: Item[]) {
  return collection.filter((item) => !collectionToDelete.includes(item));
}

function resetGenesToDelete(state: State, _: PayloadAction<null>): State {
  return {
    ...state,
    genesToDelete: [],
  };
}

function selectFilters(
  state: State,
  action: PayloadAction<{
    key: keyof State["selectedFilters"];
    options: string[];
  }>
): State {
  const { key, options } = action.payload;

  const { selectedFilters } = state;

  if (isEqual(selectedFilters[key], options)) return state;

  const newSelectedFilters = {
    ...state.selectedFilters,
    [key]: options,
  };

  return {
    ...state,
    selectedFilters: newSelectedFilters,
    expandedTissueIds: selectFiltersHandleExpandedTissueIds(state, action),
  };

  function selectFiltersHandleExpandedTissueIds(
    state: State,
    action: PayloadAction<{
      key: keyof State["selectedFilters"];
      options: string[];
    }>
  ): State["expandedTissueIds"] {
    const { expandedTissueIds, selectedFilters } = state;

    const {
      payload: { key, options: tissueIds },
    } = action;

    if (key !== "tissues") return expandedTissueIds;

    const prevTissueIds = selectedFilters.tissues;

    /**
     * (thuang): Only expand the newly added tissues, so if the user has collapsed
     * some tissues during tissue filter mode, they won't be expanded when they
     * add more tissues
     */
    const addedIds = tissueIds.filter(
      (tissueId) => !prevTissueIds.includes(tissueId)
    );

    /**
     * Entering tissue filter state: Expand all selected tissues
     */
    if (!prevTissueIds.length && tissueIds.length) {
      return tissueIds;
    }

    /**
     * Exiting tissue filter state: Collapse all tissues
     */
    if (prevTissueIds.length && !tissueIds.length) {
      return EMPTY_ARRAY;
    }

    return [
      /**
       * (thuang): Remove any expanded tissues that are no longer selected
       */
      ...expandedTissueIds.filter((expandedTissueId) =>
        tissueIds.includes(expandedTissueId)
      ),
      ...addedIds,
    ];
  }
}

function setSnapshotId(
  state: State,
  action: PayloadAction<State["snapshotId"]>
): State {
  const { payload } = action;

  return {
    ...state,
    snapshotId: payload,
  };
}

function addSelectedGenes(
  state: State,
  action: PayloadAction<State["selectedGenes"]>
): State {
  const { payload } = action;
  let newSelectedGenes = state.selectedGenes;

  newSelectedGenes = newSelectedGenes.filter((gene) => !payload.includes(gene));

  return {
    ...state,
    selectedGenes: [...payload, ...newSelectedGenes],
  };
}

export interface AddCellInfoCellTypePayload {
  cellType: CellType;
  tissueID: string;
}

function addCellInfoCellType(
  state: State,
  action: PayloadAction<AddCellInfoCellTypePayload>
): State {
  const { payload } = action;

  // Type safety, this should never happen
  if (!state.selectedOrganismId) return state;

  const newCellInfoCellType = {
    cellType: payload.cellType,
    organismID: state.selectedOrganismId,
    tissueID: payload.tissueID,
  };
  return {
    ...state,
    cellInfoCellType: newCellInfoCellType,
    geneInfoGene: null,
  };
}

function clearCellInfoCellType(state: State, _: PayloadAction<null>): State {
  return {
    ...state,
    cellInfoCellType: null,
  };
}

function clearGeneInfoGene(state: State, _: PayloadAction<null>): State {
  return {
    ...state,
    geneInfoGene: null,
  };
}

function addGeneInfoGene(state: State, action: PayloadAction<string>): State {
  const { payload } = action;

  return {
    ...state,
    geneInfoGene: payload,
  };
}

function selectGeneInfoFromXAxis(
  state: State,
  action: PayloadAction<string>
): State {
  const { payload } = action;

  return {
    ...state,
    geneInfoGene: payload,
    cellInfoCellType: null,
  };
}

function closeRightSidebar(state: State, _: PayloadAction<null>): State {
  return {
    ...state,
    geneInfoGene: null,
    cellInfoCellType: null,
  };
}

export interface LoadStateFromURLPayload {
  compare: State["compare"];
  filters: Partial<State["selectedFilters"]>;
  organism: State["selectedOrganismId"];
  tissues?: State["selectedTissues"];
  genes: State["selectedGenes"];
  cellTypes?: State["filteredCellTypes"];
}

function loadStateFromURL(
  state: State,
  action: PayloadAction<LoadStateFromURLPayload>
): State {
  const { payload } = action;

  const { compare, filters, genes, tissues } = payload;

  return {
    ...state,
    compare,
    selectedFilters: { ...state.selectedFilters, ...filters },
    selectedGenes: genes,
    selectedTissues: tissues,
    selectedOrganismId: payload.organism,
    filteredCellTypes: payload.cellTypes ?? [],
    expandedTissueIds: filters.tissues ?? [],
  };
}

function selectCompare(
  state: State,
  action: PayloadAction<State["compare"]>
): State {
  return {
    ...state,
    compare: action.payload,
  };
}

function setXAxisHeight(state: State, action: PayloadAction<number>): State {
  return {
    ...state,
    xAxisHeight: action.payload + GENE_SEARCH_BAR_HEIGHT_PX,
  };
}

export interface SetFilteredCellTypesPayload {
  filteredCellTypes: State["filteredCellTypes"];
  filteredCellTypeIds: State["filteredCellTypeIds"];
  displayedTissueIds: string[];
}

function setFilteredCellTypes(
  state: State,
  action: PayloadAction<SetFilteredCellTypesPayload>
): State {
  const {
    payload: { filteredCellTypes, filteredCellTypeIds },
  } = action;

  track(EVENTS.WMG_SELECT_CELL_TYPE, {
    cell_types: filteredCellTypeIds,
  });

  return {
    ...state,
    filteredCellTypes,
    filteredCellTypeIds,
    expandedTissueIds: setFilteredCellTypesHandleExpandedTissueIds(
      state,
      action
    ),
  };

  function setFilteredCellTypesHandleExpandedTissueIds(
    state: State,
    action: PayloadAction<SetFilteredCellTypesPayload>
  ): State["expandedTissueIds"] {
    const { filteredCellTypeIds: prevFilteredCellTypeIds } = state;

    const {
      payload: { displayedTissueIds, filteredCellTypeIds },
    } = action;

    /**
     * (thuang): Exiting cell type filter state: Collapse all tissues
     */
    if (prevFilteredCellTypeIds.length && filteredCellTypeIds.length === 0) {
      return EMPTY_ARRAY;
    }

    return displayedTissueIds;
  }
}

function toggleExpandedTissueId(
  state: State,
  action: PayloadAction<{ tissueId: string; tissueName: string }>
) {
  const {
    payload: { tissueId, tissueName },
  } = action;

  const { expandedTissueIds } = state;

  if (expandedTissueIds.includes(tissueId)) {
    return {
      ...state,
      expandedTissueIds: expandedTissueIds.filter(
        (expandedTissueId) => expandedTissueId !== tissueId
      ),
    };
  }

  track(EVENTS.WMG_TISSUE_EXPAND, { tissue: tissueName });

  return {
    ...state,
    expandedTissueIds: [...expandedTissueIds, tissueId],
  };
}

function autoExpandTissues(
  state: State,
  action: PayloadAction<string[]>
): State {
  const { payload } = action;

  return {
    ...state,
    expandedTissueIds: payload,
  };
}
