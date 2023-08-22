/**
 * Create a corresponding action creator function for each action type defined
 * in `./reducer.ts`.
 */

import { REDUCERS, State } from "./reducer";

export function setCellGuideTitle(
  payload: State["cellGuideTitle"]
): GetActionTypeOfReducer<typeof REDUCERS["setCellGuideTitle"]> {
  return {
    payload,
    type: "setCellGuideTitle",
  };
}

export function setCellGuideNav(
  payload: State["cellGuideNav"]
): GetActionTypeOfReducer<typeof REDUCERS["setCellGuideNav"]> {
  return {
    payload,
    type: "setCellGuideNav",
  };
}

type GetActionTypeOfReducer<T> = T extends (
  state: never,
  action: infer Action
) => unknown
  ? Action
  : never;
