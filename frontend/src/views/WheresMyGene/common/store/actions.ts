/**
 * Create a corresponding action creator function for each action type defined
 * in `./reducer.ts`.
 */

import { REDUCERS, State } from "./reducer";

export function deleteSelectedGenesAndSelectedCellTypeIds({
  genesToDelete,
  cellTypeIdsToDelete,
}: {
  genesToDelete: string[];
  cellTypeIdsToDelete: string[];
}): GetActionTypeOfReducer<
  typeof REDUCERS["deleteSelectedGenesAndSelectedCellTypeIds"]
> {
  return {
    payload: {
      cellTypeIdsToDelete,
      genesToDelete,
    },
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
  cellTypeIdToDelete: string
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

type GetActionTypeOfReducer<T> = T extends (
  state: never,
  action: infer Action
) => unknown
  ? Action
  : never;
