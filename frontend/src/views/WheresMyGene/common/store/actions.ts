/**
 * Create a corresponding action creator function for each action type defined
 * in `./reducer.ts`.
 */

import { CellTypeMetadata } from "../../components/HeatMap/utils";
import { CellTypeSummary } from "../types";
import { REDUCERS, State } from "./reducer";

export function deleteSelectedGenesAndSelectedCellTypeIds(): GetActionTypeOfReducer<
  typeof REDUCERS["deleteSelectedGenesAndSelectedCellTypeIds"]
> {
  return {
    payload: null,
    type: "deleteSelectedGenesAndSelectedCellTypeIds",
  };
}

export function toggleGeneToDelete(
  geneToDelete: string
): GetActionTypeOfReducer<typeof REDUCERS["toggleGeneToDelete"]> {
  return {
    payload: geneToDelete,
    type: "toggleGeneToDelete",
  };
}

export function toggleCellTypeIdToDelete(
  cellTypeIdToDelete: CellTypeMetadata
): GetActionTypeOfReducer<typeof REDUCERS["toggleCellTypeIdToDelete"]> {
  return {
    payload: cellTypeIdToDelete,
    type: "toggleCellTypeIdToDelete",
  };
}

export function selectOrganism(
  organism: State["selectedOrganism"]
): GetActionTypeOfReducer<typeof REDUCERS["selectOrganism"]> {
  return {
    payload: organism,
    type: "selectOrganism",
  };
}

export function selectGenes(
  genes: State["selectedGenes"]
): GetActionTypeOfReducer<typeof REDUCERS["selectGenes"]> {
  return {
    payload: genes,
    type: "selectGenes",
  };
}

export function selectCellTypeIds(
  cellTypeIndices: State["selectedCellTypeIds"]
): GetActionTypeOfReducer<typeof REDUCERS["selectCellTypeIds"]> {
  return {
    payload: cellTypeIndices,
    type: "selectCellTypeIds",
  };
}

export function selectTissues(
  tissues: State["selectedTissues"]
): GetActionTypeOfReducer<typeof REDUCERS["selectTissues"]> {
  return {
    payload: tissues,
    type: "selectTissues",
  };
}

export function resetGenesToDeleteAndCellTypeIdsToDelete(): GetActionTypeOfReducer<
  typeof REDUCERS["resetGenesToDeleteAndCellTypeIdsToDelete"]
> {
  return {
    payload: null,
    type: "resetGenesToDeleteAndCellTypeIdsToDelete",
  };
}

export function tissueCellTypesFetched(
  tissue: string,
  cellTypeSummaries: CellTypeSummary[]
): GetActionTypeOfReducer<typeof REDUCERS["tissueCellTypesFetched"]> {
  return {
    payload: { cellTypeSummaries, tissue },
    type: "tissueCellTypesFetched",
  };
}

export function resetTissueCellTypes(
  tissue: string,
  cellTypeSummaries: CellTypeSummary[]
): GetActionTypeOfReducer<typeof REDUCERS["resetTissueCellTypes"]> {
  return {
    payload: { cellTypeSummaries, tissue },
    type: "resetTissueCellTypes",
  };
}

type GetActionTypeOfReducer<T> = T extends (
  state: never,
  action: infer Action
) => unknown
  ? Action
  : never;
