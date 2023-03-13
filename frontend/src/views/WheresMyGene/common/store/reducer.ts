import cloneDeep from "lodash/cloneDeep";
import isEqual from "lodash/isEqual";
import {
  EMPTY_FILTERS,
  SelectedFilters,
} from "src/common/queries/wheresMyGene";
import {
  CellTypeMetadata,
  deserializeCellTypeMetadata,
} from "../../components/HeatMap/utils";
import { CellType, SORT_BY, Tissue } from "../types";
export interface PayloadAction<Payload> {
  type: keyof typeof REDUCERS;
  payload: Payload;
}
export interface State {
  cellTypeIdsToDelete: CellTypeMetadata[];
  genesToDelete: string[];
  selectedGenes: string[];
  selectedCellTypeIds: {
    [tissue: Tissue]: string[];
  };
  selectedOrganismId: string | null;
  selectedTissues: string[];
  selectedFilters: SelectedFilters;
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
}

// (thuang): If you have derived states based on the state, use `useMemo`
// to cache the derived states instead of putting them in the state.
export const INITIAL_STATE: State = {
  cellInfoCellType: null,
  cellTypeIdsToDelete: [],
  genesToDelete: [],
  selectedCellTypeIds: {},
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
};

export const REDUCERS = {
  addCellInfoCellType,
  addSelectedGenes,
  deleteSelectedGenesAndSelectedCellTypeIds,
  deleteSingleGene,
  resetGenesToDeleteAndCellTypeIdsToDelete,
  resetTissueCellTypes,
  selectCellTypeIds,
  selectFilters,
  selectGenes,
  selectOrganism,
  selectSortBy,
  selectTissues,
  setSnapshotId,
  tissueCellTypesFetched,
  toggleCellTypeIdToDelete,
  toggleGeneToDelete,
  loadStateFromURL,
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
  };
}

function deleteSelectedGenesAndSelectedCellTypeIds(
  state: State,
  _: PayloadAction<null>
): State {
  const { genesToDelete, cellTypeIdsToDelete } = state;

  if (!genesToDelete.length && !cellTypeIdsToDelete.length) {
    return state;
  }

  const { selectedGenes, selectedCellTypeIds } = state;

  const newSelectedGenes = genesToDelete.length
    ? deleteByItems<State["selectedGenes"][number]>(
        selectedGenes,
        genesToDelete
      )
    : selectedGenes;

  const newSelectedCellTypeIds = cellTypeIdsToDelete.length
    ? deleteSelectedCellTypeIdsByMetadata(
        selectedCellTypeIds,
        cellTypeIdsToDelete
      )
    : selectedCellTypeIds;

  return {
    ...state,
    cellTypeIdsToDelete: [],
    genesToDelete: [],
    selectedCellTypeIds: newSelectedCellTypeIds,
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
    selectedTissues: [],
  };
}

function selectGenes(
  state: State,
  action: PayloadAction<State["selectedGenes"]>
): State {
  return {
    ...state,
    cellTypeIdsToDelete: [],
    genesToDelete: [],
    selectedGenes: action.payload,
  };
}

function selectCellTypeIds(
  state: State,
  action: PayloadAction<State["selectedCellTypeIds"]>
): State {
  return {
    ...state,
    selectedCellTypeIds: action.payload,
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

function toggleCellTypeIdToDelete(
  state: State,
  action: PayloadAction<CellTypeMetadata>
): State {
  if (state.cellTypeIdsToDelete.includes(action.payload)) {
    return {
      ...state,
      cellTypeIdsToDelete: deleteByItems<CellTypeMetadata>(
        state.cellTypeIdsToDelete,
        [action.payload]
      ),
    };
  }

  return {
    ...state,
    cellTypeIdsToDelete: [...state.cellTypeIdsToDelete, action.payload],
  };
}

function deleteByItems<Item>(collection: Item[], collectionToDelete: Item[]) {
  return collection.filter((item) => !collectionToDelete.includes(item));
}

function resetGenesToDeleteAndCellTypeIdsToDelete(
  state: State,
  _: PayloadAction<null>
): State {
  return {
    ...state,
    cellTypeIdsToDelete: [],
    genesToDelete: [],
  };
}

function deleteSelectedCellTypeIdsByMetadata(
  selectedCellTypeIds: State["selectedCellTypeIds"],
  cellTypeMetadata: CellTypeMetadata[]
): State["selectedCellTypeIds"] {
  const newSelectedCellTypeIds = cloneDeep(selectedCellTypeIds);

  const cellTypeIdsToDeleteByTissue = cellTypeMetadata.reduce(
    (memo, metadata) => {
      const { tissue, id } = deserializeCellTypeMetadata(metadata);
      const cellTypeIds = memo[tissue] || [];
      memo[tissue] = [...cellTypeIds, id];

      return memo;
    },
    {} as State["selectedCellTypeIds"]
  );

  for (const [tissue, cellTypeIdsToDelete] of Object.entries(
    cellTypeIdsToDeleteByTissue
  )) {
    const tissueCellTypeIds = newSelectedCellTypeIds[tissue] || [];

    newSelectedCellTypeIds[tissue] = tissueCellTypeIds.filter(
      (id) => !cellTypeIdsToDelete.includes(id)
    );
  }

  return newSelectedCellTypeIds;
}

function tissueCellTypesFetched(
  state: State,
  action: PayloadAction<{
    tissue: Tissue;
    cellTypes: CellType[];
  }>
): State {
  const { tissue, cellTypes } = action.payload;

  const newCellTypeIds = cellTypes.map((cellType) => cellType.id);

  const { selectedCellTypeIds } = state;

  const oldCellTypeIds = selectedCellTypeIds[tissue];

  return {
    ...state,
    selectedCellTypeIds: {
      ...selectedCellTypeIds,
      [tissue]: oldCellTypeIds || newCellTypeIds,
    },
  };
}

function resetTissueCellTypes(
  state: State,
  action: PayloadAction<{
    tissue: Tissue;
    cellTypes: CellType[];
  }>
): State {
  const { tissue, cellTypes } = action.payload;

  const newCellTypeIds = cellTypes.map((cellType) => cellType.id);

  const { selectedCellTypeIds } = state;

  return {
    ...state,
    selectedCellTypeIds: {
      ...selectedCellTypeIds,
      [tissue]: newCellTypeIds,
    },
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
  };
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
  };
}

export interface LoadStateFromURLPayload {
  filters: Partial<State["selectedFilters"]>;
  tissues: State["selectedTissues"];
  genes: State["selectedGenes"];
}

function loadStateFromURL(
  state: State,
  action: PayloadAction<LoadStateFromURLPayload>
): State {
  const { payload } = action;

  return {
    ...state,
    selectedFilters: { ...state.selectedFilters, ...payload.filters },
    selectedTissues: payload.tissues,
    selectedGenes: payload.genes,
  };
}
