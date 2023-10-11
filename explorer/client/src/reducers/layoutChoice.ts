/*
we have a UI heuristic to pick the default layout, based on assumptions
about commonly used names.  Preferentially, pick in the following order:

  1. "umap"
  2. "tsne"
  3. "pca"
  4. give up, use the first available
*/

import type { Action, AnyAction } from "redux";
import type { RootState } from ".";
import { EmbeddingSchema, Schema } from "../common/types/schema";

function bestDefaultLayout(layouts: Array<string>): string {
  const preferredNames = ["umap", "tsne", "pca"];
  const idx = preferredNames.findIndex((name) => layouts.indexOf(name) !== -1);
  if (idx !== -1) return preferredNames[idx];
  return layouts[0];
}

function setToDefaultLayout(schema: Schema): LayoutChoiceState {
  const available = schema.layout.obs
    .map((v: EmbeddingSchema) => v.name)
    .sort();
  const current = bestDefaultLayout(available);
  const currentDimNames = schema.layout.obsByName[current].dims;
  return { available, current, currentDimNames };
}

export interface LayoutChoiceState {
  available: Array<string>;
  current: string;
  currentDimNames: Array<string>;
}

export interface LayoutChoiceAction extends Action<string> {
  layoutChoice: string;
}

const LayoutChoice = (
  state: LayoutChoiceState,
  action: AnyAction,
  nextSharedState: RootState,
): LayoutChoiceState => {
  switch (action.type) {
    case "initial data load complete": {
      // set default to default
      const { annoMatrix } = nextSharedState;
      return {
        ...state,
        ...setToDefaultLayout(annoMatrix.schema),
      };
    }

    case "set layout choice": {
      const { schema } = nextSharedState.annoMatrix;
      const current = (action as LayoutChoiceAction).layoutChoice;
      const currentDimNames = schema.layout.obsByName[current].dims;
      return { ...state, current, currentDimNames };
    }

    default: {
      return state;
    }
  }
};

export default LayoutChoice;
