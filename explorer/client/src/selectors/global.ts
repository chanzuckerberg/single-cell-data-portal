/* App dependencies */
import { selectIsUserStateDirty as annoMatrix } from "./annoMatrix";
import { selectIsUserStateDirty as categoricalSelection } from "./categoricalSelection";
import { selectIsUserStateDirty as colors } from "./colors";
import { selectIsUserStateDirty as continuousSelection } from "./continuousSelection";
import { selectIsUserStateDirty as controls } from "./controls";
import { selectIsUserStateDirty as genesets } from "./genesets";
import { RootState } from "../reducers";

/*
 Returns true if app controls have been touched indicating work is in progress.
 Touched controls include:
 - user defined category,
 - categorical selection,
 - categorical metadata color selection,
 - histogram brush,
 - individual genes selection, or
 - gene set selection.
 @params state
 @returns boolean
 */
export const selectIsUserStateDirty = (state: RootState): boolean => {
  const selectors = [
    annoMatrix,
    categoricalSelection,
    colors,
    continuousSelection,
    controls,
    genesets,
  ];

  return selectors.some((selector) => selector(state));
};
