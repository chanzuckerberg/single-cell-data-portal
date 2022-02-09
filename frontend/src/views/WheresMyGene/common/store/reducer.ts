export interface PayloadAction<Payload> {
  type: keyof typeof REDUCERS;
  payload: Payload;
}

export interface State {
  cellTypeIdsToDelete: string[];
  genesToDelete: string[];
  selectedGenes: string[];
  selectedCellTypeIds: string[];
  selectedOrganism: string;
  selectedTissues: string[];
}

export const INITIAL_STATE: State = {
  cellTypeIdsToDelete: [],
  genesToDelete: [],
  selectedCellTypeIds: [],
  selectedGenes: [],
  selectedOrganism: "",
  selectedTissues: [],
};

export const REDUCERS = {
  deleteSelectedGenesAndSelectedCellTypeIds,
  resetGenesToDeleteAndCellTypeIdsToDelete,
  selectCellTypeIds,
  selectGenes,
  selectOrganism,
  selectTissues,
  toggleCellTypeIdToDelete,
  toggleGeneToDelete,
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

function deleteSelectedGenesAndSelectedCellTypeIds(
  state: State,
  action: PayloadAction<{
    genesToDelete: string[];
    cellTypeIdsToDelete: string[];
  }>
): State {
  const { genesToDelete, cellTypeIdsToDelete } = action.payload;

  if (!genesToDelete.length && !cellTypeIdsToDelete.length) {
    return state;
  }

  const { selectedGenes, selectedCellTypeIds } = state;

  const newSelectedGenes = deleteByIds<State["selectedGenes"][number]>(
    selectedGenes,
    genesToDelete
  );

  const newSelectedCellTypeIds = deleteByIds<
    State["selectedCellTypeIds"][number]
  >(selectedCellTypeIds, cellTypeIdsToDelete);

  return {
    ...state,
    selectedCellTypeIds: newSelectedCellTypeIds,
    selectedGenes: newSelectedGenes,
  };
}

function selectOrganism(state: State, action: PayloadAction<string>): State {
  return {
    ...state,
    selectedOrganism: action.payload,
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

function toggleGeneToDelete(
  state: State,
  action: PayloadAction<string>
): State {
  if (state.genesToDelete.includes(action.payload)) {
    return {
      ...state,
      genesToDelete: deleteByIds<string>(state.genesToDelete, [action.payload]),
    };
  }

  return {
    ...state,
    genesToDelete: [...state.genesToDelete, action.payload],
  };
}

function toggleCellTypeIdToDelete(
  state: State,
  action: PayloadAction<string>
): State {
  if (state.cellTypeIdsToDelete.includes(action.payload)) {
    return {
      ...state,
      cellTypeIdsToDelete: deleteByIds<string>(state.cellTypeIdsToDelete, [
        action.payload,
      ]),
    };
  }

  return {
    ...state,
    cellTypeIdsToDelete: [...state.cellTypeIdsToDelete, action.payload],
  };
}

function deleteByIds<Item>(collection: Item[], collectionToDelete: Item[]) {
  return collection.filter((item) => !collectionToDelete.includes(item));
}

function resetGenesToDeleteAndCellTypeIdsToDelete(
  state: State,
  // eslint-disable-next-line @typescript-eslint/no-unused-vars
  _: PayloadAction<null>
): State {
  return {
    ...state,
    cellTypeIdsToDelete: [],
    genesToDelete: [],
  };
}
