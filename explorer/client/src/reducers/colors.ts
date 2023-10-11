/*
 * Color By UI state
 */

import { Action, AnyAction } from "redux";

export interface UserColor {
  colors: { [label: string]: [number, number, number] };
  scale: (label: string) => d3.RGBColor;
}

export interface ConvertedUserColors {
  [category: string]: UserColor;
}

export interface ColorsState {
  // by continuous, by expression
  colorMode: Action["type"] | null;
  // tissue, Apod
  colorAccessor: string | null;
  userColors?: ConvertedUserColors;
}

const ColorsReducer = (
  state: ColorsState = {
    colorMode: null,
    colorAccessor: null,
  },
  action: AnyAction,
): ColorsState => {
  switch (action.type) {
    case "universe: user color load success": {
      const { userColors } = action;
      return {
        ...state,
        userColors,
      };
    }

    case "annotation: category edited": {
      const { colorAccessor } = state;
      if (action.metadataField !== colorAccessor) {
        return state;
      }
      /* else update colorAccessor */
      return {
        ...state,
        colorAccessor: action.newCategoryText,
      };
    }

    case "annotation: delete category": {
      const { colorAccessor } = state;
      if (action.metadataField !== colorAccessor) {
        return state;
      }
      /* else reset */
      return {
        ...state,
        colorMode: null,
        colorAccessor: null,
      };
    }

    case "reset colorscale": {
      return {
        ...state,
        colorMode: null,
        colorAccessor: null,
      };
    }

    case "color by categorical metadata":
    case "color by continuous metadata": {
      /* toggle between this mode and reset */
      const resetCurrent =
        action.type === state.colorMode &&
        action.colorAccessor === state.colorAccessor;
      const colorMode = !resetCurrent ? action.type : null;
      const colorAccessor = !resetCurrent ? action.colorAccessor : null;

      return {
        ...state,
        colorMode,
        colorAccessor,
      };
    }

    case "color by expression": {
      /* toggle between this mode and reset */
      const resetCurrent =
        action.type === state.colorMode && action.gene === state.colorAccessor;
      const colorMode = !resetCurrent ? action.type : null;
      const colorAccessor = !resetCurrent ? action.gene : null;

      return {
        ...state,
        colorMode,
        colorAccessor,
      };
    }

    case "color by geneset mean expression": {
      /* toggle between this mode and reset */
      const resetCurrent =
        action.type === state.colorMode &&
        action.geneset === state.colorAccessor;
      const colorMode = !resetCurrent ? action.type : null;
      const colorAccessor = !resetCurrent ? action.geneset : null;

      return {
        ...state,
        colorMode,
        colorAccessor,
      };
    }

    default: {
      return state;
    }
  }
};

export default ColorsReducer;
