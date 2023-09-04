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
  cellTypes: string[];
}

export type QueryGroup = Filters;

export interface QueryGroups {
  queryGroup1: QueryGroup;
  queryGroup2: QueryGroup;
}

export type QueryGroupsWithNames = QueryGroups;
export interface State {
  organismId: string | null;
  queryGroups: QueryGroups;
  queryGroupsWithNames: QueryGroupsWithNames;
  submittedQueryGroups: QueryGroups | null;
  submittedQueryGroupsWithNames: QueryGroupsWithNames | null;
  snapshotId: string | null;
}

export const EMPTY_FILTERS = {
  datasets: [],
  developmentStages: [],
  diseases: [],
  ethnicities: [],
  sexes: [],
  tissues: [],
  cellTypes: [],
};

export const INITIAL_STATE: State = {
  organismId: null,
  snapshotId: null,
  queryGroups: { queryGroup1: EMPTY_FILTERS, queryGroup2: EMPTY_FILTERS },
  queryGroupsWithNames: {
    queryGroup1: EMPTY_FILTERS,
    queryGroup2: EMPTY_FILTERS,
  },
  submittedQueryGroups: null,
  submittedQueryGroupsWithNames: null,
};

export const REDUCERS = {
  selectOrganism,
  setSnapshotId,
  selectQueryGroup1Filters,
  selectQueryGroup2Filters,
  setQueryGroup1Filters,
  setQueryGroup2Filters,
  clearQueryGroup1Filters,
  clearQueryGroup2Filters,
  copyCellGroup1,
  submitQueryGroups,
  clearSubmittedQueryGroups,
};

function setQueryGroup1Filters(
  state: State,
  action: PayloadAction<QueryGroup>
): State {
  const { payload } = action;

  return {
    ...state,
    queryGroups: { ...state.queryGroups, queryGroup1: payload },
  };
}

function setQueryGroup2Filters(
  state: State,
  action: PayloadAction<QueryGroup>
): State {
  const { payload } = action;

  return {
    ...state,
    queryGroups: { ...state.queryGroups, queryGroup2: payload },
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

function selectQueryGroup1Filters(
  state: State,
  action: PayloadAction<{
    key: keyof QueryGroup;
    options: { id: string; name: string }[];
  }>
): State {
  const { key, options } = action.payload;

  const { queryGroups, queryGroupsWithNames } = state;

  const { queryGroup1 } = queryGroups;

  const newQueryGroup = { ...queryGroup1 };
  newQueryGroup[key] = options.map((option) => option.id);

  const { queryGroup1: queryGroupWithNames1 } = queryGroupsWithNames;

  const newQueryGroupWithNames = { ...queryGroupWithNames1 };
  newQueryGroupWithNames[key] = options.map((option) => option.name);
  return {
    ...state,
    queryGroups: { ...queryGroups, queryGroup1: newQueryGroup },
    queryGroupsWithNames: {
      ...queryGroupsWithNames,
      queryGroup1: newQueryGroupWithNames,
    },
  };
}

function selectQueryGroup2Filters(
  state: State,
  action: PayloadAction<{
    key: keyof QueryGroup;
    options: { id: string; name: string }[];
  }>
): State {
  const { key, options } = action.payload;

  const { queryGroups, queryGroupsWithNames } = state;

  const { queryGroup2 } = queryGroups;

  const newQueryGroup = { ...queryGroup2 };
  newQueryGroup[key] = options.map((option) => option.id);

  const { queryGroup2: queryGroupWithNames2 } = queryGroupsWithNames;

  const newQueryGroupWithNames = { ...queryGroupWithNames2 };
  newQueryGroupWithNames[key] = options.map((option) => option.name);

  return {
    ...state,
    queryGroups: { ...queryGroups, queryGroup2: newQueryGroup },
    queryGroupsWithNames: {
      ...queryGroupsWithNames,
      queryGroup2: newQueryGroupWithNames,
    },
  };
}

function submitQueryGroups(state: State, _: PayloadAction<null>): State {
  const { queryGroups, queryGroupsWithNames } = state;

  return {
    ...state,
    submittedQueryGroups: queryGroups,
    submittedQueryGroupsWithNames: queryGroupsWithNames,
  };
}

function clearSubmittedQueryGroups(
  state: State,
  _: PayloadAction<null>
): State {
  return {
    ...state,
    submittedQueryGroups: null,
    submittedQueryGroupsWithNames: null,
  };
}

function clearQueryGroup1Filters(state: State, _: PayloadAction<null>): State {
  const { queryGroups, queryGroupsWithNames } = state;

  return {
    ...state,
    queryGroups: { ...queryGroups, queryGroup1: EMPTY_FILTERS },
    queryGroupsWithNames: {
      ...queryGroupsWithNames,
      queryGroup1: EMPTY_FILTERS,
    },
  };
}

function clearQueryGroup2Filters(state: State, _: PayloadAction<null>): State {
  const { queryGroups, queryGroupsWithNames } = state;

  return {
    ...state,
    queryGroups: { ...queryGroups, queryGroup2: EMPTY_FILTERS },
    queryGroupsWithNames: {
      ...queryGroupsWithNames,
      queryGroup2: EMPTY_FILTERS,
    },
  };
}

function copyCellGroup1(state: State, _: PayloadAction<null>): State {
  const { queryGroups, queryGroupsWithNames } = state;

  return {
    ...state,
    queryGroups: { ...queryGroups, queryGroup2: queryGroups.queryGroup1 },
    queryGroupsWithNames: {
      ...queryGroupsWithNames,
      queryGroup2: queryGroupsWithNames.queryGroup1,
    },
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
