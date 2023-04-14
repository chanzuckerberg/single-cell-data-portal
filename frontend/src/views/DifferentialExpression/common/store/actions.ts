import { QueryGroup, REDUCERS, State } from "./reducer";

export function selectOrganism(
  organismId: State["organismId"]
): GetActionTypeOfReducer<typeof REDUCERS["selectOrganism"]> {
  return {
    payload: organismId,
    type: "selectOrganism",
  };
}

export function selectFilters(
  key: keyof State["selectedFilters"],
  options: string[]
): GetActionTypeOfReducer<typeof REDUCERS["selectFilters"]> {
  return {
    payload: { key, options },
    type: "selectFilters",
  };
}

export function setSelectedFilterNames(
  key: keyof State["selectedFilterNames"],
  options: string[]
): GetActionTypeOfReducer<typeof REDUCERS["setSelectedFilterNames"]> {
  return {
    payload: { key, options },
    type: "setSelectedFilterNames",
  };
}

export function setSnapshotId(
  snapshotId: State["snapshotId"]
): GetActionTypeOfReducer<typeof REDUCERS["setSnapshotId"]> {
  return {
    payload: snapshotId,
    type: "setSnapshotId",
  };
}

export function addQueryGroup(
  queryGroup: QueryGroup | null
): GetActionTypeOfReducer<typeof REDUCERS["addQueryGroup"]> {
  return {
    payload: queryGroup,
    type: "addQueryGroup",
  };
}

type GetActionTypeOfReducer<T> = T extends (
  state: never,
  action: infer Action
) => unknown
  ? Action
  : never;
