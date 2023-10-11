import { AnyAction } from "redux";
import type { RootState } from ".";
import { ControlsHelpers as CH } from "../util/stateManager";

/*
State is an object, with a key for each categorical annotation, and a
value which is a Map of label->t/f, reflecting selection state for the label.

Label state default (if missing) is up to the component, but typically true.

{
  "louvain": Map(),
  ...
}
*/
interface CategoricalSelectionState {
  [key: string]: Map<string, boolean>;
}
const CategoricalSelection = (
  state: CategoricalSelectionState,
  action: AnyAction,
  nextSharedState: RootState,
): CategoricalSelectionState => {
  switch (action.type) {
    case "initial data load complete":
    case "subset to selection":
    case "reset subset":
    case "set clip quantiles": {
      const { annoMatrix } = nextSharedState;
      const newState = CH.createCategoricalSelection(
        CH.selectableCategoryNames(annoMatrix.schema),
      );
      return newState;
    }

    case "categorical metadata filter select":
    case "categorical metadata filter deselect":
    case "categorical metadata filter none of these":
    case "categorical metadata filter all of these": {
      const { metadataField, labelSelectionState } = action;
      return {
        ...state,
        [metadataField]: labelSelectionState,
      };
    }

    default: {
      return state;
    }
  }
};

export default CategoricalSelection;
