/*
Reducer for the obsCrossfilter
*/

import { AnyAction } from "redux";
import { AnnoMatrixObsCrossfilter } from "../annoMatrix";

const ObsCrossfilter = (
  state: AnnoMatrixObsCrossfilter | null = null,
  action: AnyAction,
) => {
  if (action.obsCrossfilter) {
    return action.obsCrossfilter;
  }
  return state;
};

export default ObsCrossfilter;
