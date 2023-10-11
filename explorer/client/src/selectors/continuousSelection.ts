/* App dependencies */
import { RootState } from "../reducers";
import { ContinuousSelectionState } from "../reducers/continuousSelection";

export const selectContinuousSelection = (
  state: RootState,
): ContinuousSelectionState => state.continuousSelection;

/*
 Returns true if histogram brush controls have been touched indicating work is in progress.
 @param continuousSelection from state
 @returns boolean
 */
export const selectIsUserStateDirty = (
  state: ContinuousSelectionState,
): boolean => {
  const continuousSelection = selectContinuousSelection(state);

  return Boolean(Object.keys(continuousSelection).length);
};
