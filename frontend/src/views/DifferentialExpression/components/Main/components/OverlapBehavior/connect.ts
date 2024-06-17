import { useContext } from "react";
import {
  DispatchContext,
  StateContext,
} from "src/views/DifferentialExpression/common/store";
import { setExcludeOverlappingCells } from "src/views/DifferentialExpression/common/store/actions";
import { ExcludeOverlappingCells } from "src/views/DifferentialExpression/common/types";

export const useConnect = () => {
  const { excludeOverlappingCells: activeValue } = useContext(StateContext);
  const dispatch = useContext(DispatchContext);

  const setActiveValue = (excludeOverlappingCells: ExcludeOverlappingCells) => {
    if (!dispatch) return;
    dispatch(setExcludeOverlappingCells(excludeOverlappingCells));
  };
  return {
    activeValue,
    setActiveValue,
  };
};
