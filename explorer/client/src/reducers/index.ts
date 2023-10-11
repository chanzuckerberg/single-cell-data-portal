import {
  createStore,
  applyMiddleware,
  AnyAction,
  Action,
  Reducer,
} from "redux";
import thunk, { ThunkDispatch } from "redux-thunk";

import cascadeReducers from "./cascade";
import undoable from "./undoable";
import datasetMetadata from "./datasetMetadata";
import config from "./config";
import annoMatrix from "./annoMatrix";
import obsCrossfilter from "./obsCrossfilter";
import categoricalSelection from "./categoricalSelection";
import continuousSelection from "./continuousSelection";
import graphSelection from "./graphSelection";
import colors from "./colors";
import differential from "./differential";
import layoutChoice from "./layoutChoice";
import controls from "./controls";
import annotations from "./annotations";
import genesets from "./genesets";
import genesetsUI from "./genesetsUI";
import centroidLabels from "./centroidLabels";
import pointDialation from "./pointDilation";
import quickGenes from "./quickGenes";

import { gcMiddleware as annoMatrixGC } from "../annoMatrix";

import undoableConfig from "./undoableConfig";

const AppReducer = undoable(
  cascadeReducers([
    ["config", config],
    ["annoMatrix", annoMatrix],
    ["obsCrossfilter", obsCrossfilter],
    ["annotations", annotations],
    ["genesets", genesets],
    ["genesetsUI", genesetsUI],
    ["layoutChoice", layoutChoice],
    ["categoricalSelection", categoricalSelection],
    ["continuousSelection", continuousSelection],
    ["graphSelection", graphSelection],
    ["colors", colors],
    ["controls", controls],
    ["quickGenes", quickGenes],
    ["differential", differential],
    ["centroidLabels", centroidLabels],
    ["pointDilation", pointDialation],
    ["datasetMetadata", datasetMetadata],
  ]),
  [
    "annoMatrix",
    "obsCrossfilter",
    "categoricalSelection",
    "continuousSelection",
    "graphSelection",
    "colors",
    "controls",
    "quickGenes",
    "differential",
    "layoutChoice",
    "centroidLabels",
    "genesets",
    "annotations",
  ],
  undoableConfig,
);

const RootReducer: Reducer = (state: RootState, action: Action) => {
  // when a logout action is dispatched it will reset redux state
  if (action.type === "reset") {
    state = undefined;
  }

  return AppReducer(state, action);
};

const store = createStore(RootReducer, applyMiddleware(thunk, annoMatrixGC));

export type RootState = ReturnType<typeof store.getState>;

export type AppDispatch = ThunkDispatch<RootState, never, AnyAction>;

export type GetState = () => RootState;

export default store;
