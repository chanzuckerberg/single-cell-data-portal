/*
Reducer for the annoMatrix
*/

import { AnyAction } from "redux";
import AnnoMatrix from "../annoMatrix/annoMatrix";

const AnnoMatrixReducer = (
  state: AnnoMatrix | null = null,
  action: AnyAction,
) => {
  if (action.annoMatrix) {
    return action.annoMatrix;
  }
  return state;
};

export default AnnoMatrixReducer;
