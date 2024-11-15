import { FilterOption, QueryGroup, REDUCERS, State } from "./reducer";

export function selectOrganism(
  organismId: State["organismId"]
): GetActionTypeOfReducer<(typeof REDUCERS)["selectOrganism"]> {
  return {
    payload: organismId,
    type: "selectOrganism",
  };
}

export function setSelectedOptionsGroup1(
  key: keyof QueryGroup,
  options: FilterOption[]
): GetActionTypeOfReducer<(typeof REDUCERS)["setSelectedOptionsGroup1"]> {
  return {
    payload: { key, options },
    type: "setSelectedOptionsGroup1",
  };
}

export function setSelectedOptionsGroup2(
  key: keyof QueryGroup,
  options: FilterOption[]
): GetActionTypeOfReducer<(typeof REDUCERS)["setSelectedOptionsGroup2"]> {
  return {
    payload: { key, options },
    type: "setSelectedOptionsGroup2",
  };
}

export function selectQueryGroup1Filters(
  key: keyof QueryGroup,
  options: { id: string; name: string }[]
): GetActionTypeOfReducer<(typeof REDUCERS)["selectQueryGroup1Filters"]> {
  return {
    payload: { key, options },
    type: "selectQueryGroup1Filters",
  };
}

export function selectQueryGroup2Filters(
  key: keyof QueryGroup,
  options: { id: string; name: string }[]
): GetActionTypeOfReducer<(typeof REDUCERS)["selectQueryGroup2Filters"]> {
  return {
    payload: { key, options },
    type: "selectQueryGroup2Filters",
  };
}

export function submitQueryGroups(): GetActionTypeOfReducer<
  (typeof REDUCERS)["submitQueryGroups"]
> {
  return {
    payload: null,
    type: "submitQueryGroups",
  };
}
export function clearSubmittedQueryGroups(): GetActionTypeOfReducer<
  (typeof REDUCERS)["clearSubmittedQueryGroups"]
> {
  return {
    payload: null,
    type: "clearSubmittedQueryGroups",
  };
}

export function clearQueryGroup1Filters(): GetActionTypeOfReducer<
  (typeof REDUCERS)["clearQueryGroup1Filters"]
> {
  return {
    payload: null,
    type: "clearQueryGroup1Filters",
  };
}

export function clearQueryGroup2Filters(): GetActionTypeOfReducer<
  (typeof REDUCERS)["clearQueryGroup2Filters"]
> {
  return {
    payload: null,
    type: "clearQueryGroup2Filters",
  };
}

export function setExcludeOverlappingCells(
  excludeOverlappingCells: State["excludeOverlappingCells"]
): GetActionTypeOfReducer<(typeof REDUCERS)["setExcludeOverlappingCells"]> {
  return {
    payload: excludeOverlappingCells,
    type: "setExcludeOverlappingCells",
  };
}

export function setSnapshotId(
  snapshotId: State["snapshotId"]
): GetActionTypeOfReducer<(typeof REDUCERS)["setSnapshotId"]> {
  return {
    payload: snapshotId,
    type: "setSnapshotId",
  };
}

type GetActionTypeOfReducer<T> = T extends (
  state: never,
  action: infer Action
) => unknown
  ? Action
  : never;
