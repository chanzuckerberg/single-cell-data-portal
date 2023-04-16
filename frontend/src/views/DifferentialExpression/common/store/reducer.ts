import isEqual from "lodash/isEqual";

export interface PayloadAction<Payload> {
  type: keyof typeof REDUCERS;
  payload: Payload;
}

export interface Filters {
  datasets: string[];
  developmentStages: string[];
  diseases: string[];
  ethnicities: string[];
  sexes: string[];
  tissues: string[];
}
type FilterNames = Filters;

export interface QueryGroup extends Filters {
  cellTypes: string[];
}

export type QueryGroupWithNames = QueryGroup;
export interface State {
  organismId: string | null;
  selectedFilters: Filters;
  selectedFilterNames: FilterNames;
  queryGroups: QueryGroup[] | null;
  queryGroupsWithNames: QueryGroupWithNames[] | null;
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

const EMPTY_QUERY_GROUP = {
  cellTypes: [],
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
  queryGroups: null,
  queryGroupsWithNames: null,
};

export const REDUCERS = {
  selectOrganism,
  selectFilters,
  setSnapshotId,
  setSelectedFilterNames,
  addQueryGroup,
  selectQueryGroupFilters,
  deleteQueryGroup,
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

function addQueryGroup(state: State, _: PayloadAction<null>): State {
  const { queryGroups, queryGroupsWithNames } = state;

  const newQueryGroups = queryGroups ? Array.from(queryGroups) : [];
  newQueryGroups.push(EMPTY_QUERY_GROUP);

  const newQueryGroupsWithNames = queryGroupsWithNames
    ? Array.from(queryGroupsWithNames)
    : [];
  newQueryGroupsWithNames.push(EMPTY_QUERY_GROUP);

  return {
    ...state,
    queryGroups: newQueryGroups,
    queryGroupsWithNames: newQueryGroupsWithNames,
  };
}

function deleteQueryGroup(state: State, action: PayloadAction<number>): State {
  const { queryGroups, queryGroupsWithNames } = state;

  const newQueryGroups: QueryGroup[] = [];
  queryGroups?.forEach((queryGroup, index) => {
    if (index !== action.payload) {
      newQueryGroups.push(queryGroup);
    }
  });

  const newQueryGroupsWithNames: QueryGroupWithNames[] = [];
  queryGroupsWithNames?.forEach((queryGroupWithNames, index) => {
    if (index !== action.payload) {
      newQueryGroupsWithNames.push(queryGroupWithNames);
    }
  });

  return {
    ...state,
    queryGroups: newQueryGroups,
    queryGroupsWithNames: newQueryGroupsWithNames,
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

function selectQueryGroupFilters(
  state: State,
  action: PayloadAction<{
    key: keyof QueryGroup;
    options: { id: string; name: string }[];
    index: number;
  }>
): State {
  const { key, options, index } = action.payload;

  const { queryGroups, queryGroupsWithNames } = state;

  const newQueryGroups = queryGroups ? Array.from(queryGroups) : [];

  const newQueryGroup = { ...newQueryGroups[index] };
  newQueryGroup[key] = options.map((option) => option.id);
  newQueryGroups[index] = newQueryGroup;

  const newQueryGroupsWithNames = queryGroupsWithNames
    ? Array.from(queryGroupsWithNames)
    : [];
  const newQueryGroupWithNames = { ...newQueryGroupsWithNames[index] };
  newQueryGroupWithNames[key] = options.map((option) => option.name);
  newQueryGroupsWithNames[index] = newQueryGroupWithNames;

  return {
    ...state,
    queryGroups: newQueryGroups,
    queryGroupsWithNames: newQueryGroupsWithNames,
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
