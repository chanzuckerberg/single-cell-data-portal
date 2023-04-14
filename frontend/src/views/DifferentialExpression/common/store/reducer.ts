import isEqual from "lodash/isEqual";

export interface PayloadAction<Payload> {
  type: keyof typeof REDUCERS;
  payload: Payload;
}
export interface State {
  organismId: string | null;
  selectedFilters: {
    datasets: string[];
    developmentStages: string[];
    diseases: string[];
    ethnicities: string[];
    sexes: string[];
    tissues: string[];
  };
  selectedFilterNames: {
    datasets: string[];
    developmentStages: string[];
    diseases: string[];
    ethnicities: string[];
    sexes: string[];
    tissues: string[];
  };
  snapshotId: string | null;
}

const EMPTY_FILTERS = {
  datasets: [],
  developmentStages: [],
  diseases: [],
  ethnicities: [],
  sexes: [],
  tissues: [],
};

export const INITIAL_STATE: State = {
  organismId: null,
  selectedFilters: EMPTY_FILTERS,
  selectedFilterNames: EMPTY_FILTERS,
  snapshotId: null,
};

export const REDUCERS = {
  selectOrganism,
  selectFilters,
  setSnapshotId,
  setSelectedFilterNames,
};

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

function selectOrganism(
  state: State,
  action: PayloadAction<string | null>
): State {
  if (state.organismId === action.payload) {
    return state;
  }

  return {
    ...state,
    organismId: action.payload,
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

function setSelectedFilterNames(
  state: State,
  action: PayloadAction<{
    key: keyof State["selectedFilterNames"];
    options: string[];
  }>
): State {
  const { key, options } = action.payload;

  const { selectedFilterNames } = state;

  if (isEqual(selectedFilterNames[key], options)) return state;

  const newSelectedFilterNames = {
    ...state.selectedFilterNames,
    [key]: options,
  };

  return {
    ...state,
    selectedFilterNames: newSelectedFilterNames,
  };
}

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
