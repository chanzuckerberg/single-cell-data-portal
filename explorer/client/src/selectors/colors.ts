/* App dependencies */
import { RootState } from "../reducers";
import { ColorsState } from "../reducers/colors";

export const selectColors = (state: RootState): ColorsState => state.colors;

/*
 Returns true if categorical metadata color selection controls have been touched indicating work is in progress.
 @param colors from state
 @returns boolean
 */
export const selectIsUserStateDirty = (state: ColorsState): boolean => {
  const colors = selectColors(state);

  return Boolean(colors.colorAccessor);
};
