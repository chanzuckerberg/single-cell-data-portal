import type { Action, AnyAction } from "redux";

import { makeContinuousDimensionName } from "../util/nameCreators";

import type { ContinuousNamespace } from "../util/nameCreators";

export interface ContinuousSelectionAction extends Action<string> {
  continuousNamespace: ContinuousNamespace;
  selection: string;
  range: [number, number];
}

export interface ContinuousSelectionState {
  [name: string]: [number, number];
}

const ContinuousSelection = (
  state: ContinuousSelectionState = {},
  action: AnyAction,
): ContinuousSelectionState => {
  switch (action.type) {
    case "reset subset":
    case "subset to selection":
    case "set clip quantiles": {
      return {};
    }
    case "continuous metadata histogram start":
    case "continuous metadata histogram brush":
    case "continuous metadata histogram end": {
      const { continuousNamespace, selection, range } =
        action as ContinuousSelectionAction;

      const name = makeContinuousDimensionName(continuousNamespace, selection);
      return {
        ...state,
        [name]: range,
      };
    }
    case "continuous metadata histogram cancel": {
      const { continuousNamespace, selection } =
        action as ContinuousSelectionAction;

      const name = makeContinuousDimensionName(continuousNamespace, selection);
      const { [name]: deletedField, ...newState } = state;
      return newState;
    }
    default: {
      return state;
    }
  }
};

export default ContinuousSelection;
