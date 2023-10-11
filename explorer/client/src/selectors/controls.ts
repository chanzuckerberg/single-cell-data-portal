/* App dependencies */
import { RootState } from "../reducers";

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- update typings once controls reducer state is typed.
export const selectControls = (state: RootState): any => state.controls;

/*
 Returns true if individual genes have been created indicating work is in progress.
 @param controls from state
 @returns boolean
 */
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- update typings once controls reducer state is typed.
export const selectIsUserStateDirty = (state: any): boolean => {
  const controls = selectControls(state);

  return Boolean(controls?.userDefinedGenes?.length);
};
