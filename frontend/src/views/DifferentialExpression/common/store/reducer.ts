import { ExcludeOverlappingCells } from "../types";

export interface PayloadAction<Payload> {
  type: keyof typeof REDUCERS;
  payload: Payload;
}

export interface FilterOption {
  name: string;
  id: string;
  unavailable?: boolean;
}

export interface QueryGroup {
  developmentStages: string[];
  diseases: string[];
  ethnicities: string[];
  sexes: string[];
  tissues: string[];
  cellTypes: string[];
  publicationCitations: string[];
}

export interface QueryGroups {
  queryGroup1: QueryGroup;
  queryGroup2: QueryGroup;
}

type QueryGroupsWithNames = QueryGroups;
export interface State {
  organismId: string | null;
  queryGroups: QueryGroups;
  queryGroupsWithNames: QueryGroupsWithNames;
  submittedQueryGroups: QueryGroups | null;
  submittedQueryGroupsWithNames: QueryGroupsWithNames | null;
  snapshotId: string | null;
  excludeOverlappingCells: ExcludeOverlappingCells;
  selectedOptionsGroup1: {
    [key in keyof QueryGroup]: FilterOption[];
  };
  selectedOptionsGroup2: {
    [key in keyof QueryGroup]: FilterOption[];
  };
}

export const EMPTY_FILTERS = {
  developmentStages: [],
  diseases: [],
  ethnicities: [],
  sexes: [],
  tissues: [],
  cellTypes: [],
  publicationCitations: [],
};

export const INITIAL_STATE: State = {
  organismId: null,
  snapshotId: null,
  excludeOverlappingCells: "retainBoth",
  queryGroups: { queryGroup1: EMPTY_FILTERS, queryGroup2: EMPTY_FILTERS },
  queryGroupsWithNames: {
    queryGroup1: EMPTY_FILTERS,
    queryGroup2: EMPTY_FILTERS,
  },
  submittedQueryGroups: null,
  submittedQueryGroupsWithNames: null,
  selectedOptionsGroup1: EMPTY_FILTERS,
  selectedOptionsGroup2: EMPTY_FILTERS,
};

export const REDUCERS = {
  selectOrganism,
  setSnapshotId,
  selectQueryGroup1Filters,
  selectQueryGroup2Filters,
  clearQueryGroup1Filters,
  clearQueryGroup2Filters,
  submitQueryGroups,
  clearSubmittedQueryGroups,
  setExcludeOverlappingCells,
  setSelectedOptionsGroup1,
  setSelectedOptionsGroup2,
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

function setExcludeOverlappingCells(
  state: State,
  action: PayloadAction<ExcludeOverlappingCells>
): State {
  return {
    ...state,
    excludeOverlappingCells: action.payload,
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

function setSelectedOptionsGroup1(
  state: State,
  action: PayloadAction<{ key: keyof QueryGroup; options: FilterOption[] }>
): State {
  const { key, options } = action.payload;

  return {
    ...state,
    selectedOptionsGroup1: { ...state.selectedOptionsGroup1, [key]: options },
  };
}

function setSelectedOptionsGroup2(
  state: State,
  action: PayloadAction<{ key: keyof QueryGroup; options: FilterOption[] }>
): State {
  const { key, options } = action.payload;

  return {
    ...state,
    selectedOptionsGroup2: { ...state.selectedOptionsGroup2, [key]: options },
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
  const {
    queryGroups,
    queryGroupsWithNames,
    selectedOptionsGroup1,
    selectedOptionsGroup2,
  } = state;

  const newQueryGroup1 = { ...queryGroups.queryGroup1 };
  const newQueryGroup2 = { ...queryGroups.queryGroup2 };
  for (const key in selectedOptionsGroup1) {
    newQueryGroup1[key as keyof QueryGroup] = selectedOptionsGroup1[
      key as keyof QueryGroup
    ]
      .filter((option: FilterOption) => !option.unavailable)
      .map((option: FilterOption) => option.id);
  }
  for (const key in selectedOptionsGroup2) {
    newQueryGroup2[key as keyof QueryGroup] = selectedOptionsGroup2[
      key as keyof QueryGroup
    ]
      .filter((option: FilterOption) => !option.unavailable)
      .map((option: FilterOption) => option.id);
  }

  const newQueryGroupWithNames1 = { ...queryGroupsWithNames.queryGroup1 };
  const newQueryGroupWithNames2 = { ...queryGroupsWithNames.queryGroup2 };

  for (const key in selectedOptionsGroup1) {
    newQueryGroupWithNames1[key as keyof QueryGroup] = selectedOptionsGroup1[
      key as keyof QueryGroup
    ]
      .filter((option: FilterOption) => !option.unavailable)
      .map((option: FilterOption) => option.name);
  }
  for (const key in selectedOptionsGroup2) {
    newQueryGroupWithNames2[key as keyof QueryGroup] = selectedOptionsGroup2[
      key as keyof QueryGroup
    ]
      .filter((option: FilterOption) => !option.unavailable)
      .map((option: FilterOption) => option.name);
  }

  return {
    ...state,
    submittedQueryGroups: {
      queryGroup1: newQueryGroup1,
      queryGroup2: newQueryGroup2,
    },
    submittedQueryGroupsWithNames: {
      queryGroup1: newQueryGroupWithNames1,
      queryGroup2: newQueryGroupWithNames2,
    },
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
