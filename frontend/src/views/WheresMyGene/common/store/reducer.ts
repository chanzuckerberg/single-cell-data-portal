import isEqual from "lodash/isEqual";
import { EMPTY_FILTERS } from "src/common/queries/wheresMyGene";
import { CellType, SORT_BY } from "../types";
export interface PayloadAction<Payload> {
  type: keyof typeof REDUCERS;
  payload: Payload;
}
export interface State {
  genesToDelete: string[];
  selectedGenes: string[];
  selectedOrganismId: string | null;
  selectedTissues: string[];
  selectedFilters: {
    datasets: string[];
    developmentStages: string[];
    diseases: string[];
    ethnicities: string[];
    sexes: string[];
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
}

// (thuang): If you have derived states based on the state, use `useMemo`
// to cache the derived states instead of putting them in the state.
export const INITIAL_STATE: State = {
  cellInfoCellType: null,
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
};

export const REDUCERS = {
  addCellInfoCellType,
  addSelectedGenes,
  deleteSelectedGenes,
  deleteSingleGene,
  resetGenesToDelete,
  selectFilters,
  selectGenes,
  selectOrganism,
  selectSortBy,
  selectTissues,
  setSnapshotId,
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
    selectedTissues: [],
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
