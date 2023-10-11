/* App dependencies */
import { RootState } from "../reducers";
import { State } from "../reducers/genesets";

export const selectGeneSets = (state: RootState): State => state.genesets;

/*
 Returns true if genesets have been created indicating work is in progress.
 @param genesets from state
 @returns boolean
 */
export const selectIsUserStateDirty = (state: State): boolean => {
  const geneSets = selectGeneSets(state);

  return Boolean(geneSets?.genesets?.size);
};
