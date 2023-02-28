/**
 * Create a corresponding action creator function for each action type defined
 * in `./reducer.ts`.
 */

import { CellTypeMetadata } from "../../components/HeatMap/utils";
import { CellType, Tissue } from "../types";
import {
  AddCellInfoCellTypePayload,
  LoadStateFromURLPayload,
  REDUCERS,
  State,
} from "./reducer";

export function deleteSingleGene(
  geneToDelete: string
): GetActionTypeOfReducer<typeof REDUCERS["deleteSingleGene"]> {
  return {
    payload: geneToDelete,
    type: "deleteSingleGene",
  };
}

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
  organismId: State["selectedOrganismId"]
): GetActionTypeOfReducer<typeof REDUCERS["selectOrganism"]> {
  return {
    payload: organismId,
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

export function addSelectedGenes(
  genes: State["selectedGenes"]
): GetActionTypeOfReducer<typeof REDUCERS["addSelectedGenes"]> {
  return {
    payload: genes,
    type: "addSelectedGenes",
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

export function selectSortBy(
  sortBy: Partial<State["sortBy"]>
): GetActionTypeOfReducer<typeof REDUCERS["selectSortBy"]> {
  return {
    payload: sortBy,
    type: "selectSortBy",
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
  tissue: Tissue,
  cellTypes: CellType[]
): GetActionTypeOfReducer<typeof REDUCERS["tissueCellTypesFetched"]> {
  return {
    payload: { cellTypes, tissue },
    type: "tissueCellTypesFetched",
  };
}

export function resetTissueCellTypes(
  tissue: Tissue,
  cellTypes: CellType[]
): GetActionTypeOfReducer<typeof REDUCERS["resetTissueCellTypes"]> {
  return {
    payload: { cellTypes, tissue },
    type: "resetTissueCellTypes",
  };
}

export function selectFilters(
  key: keyof State["selectedFilters"],
  options: string[]
): GetActionTypeOfReducer<typeof REDUCERS["selectFilters"]> {
  return {
    payload: { key, options },
    type: "selectFilters",
  };
}

export function setSnapshotId(
  snapshotId: State["snapshotId"]
): GetActionTypeOfReducer<typeof REDUCERS["setSnapshotId"]> {
  return {
    payload: snapshotId,
    type: "setSnapshotId",
  };
}

export function addCellInfoCellType(
  payload: AddCellInfoCellTypePayload
): GetActionTypeOfReducer<typeof REDUCERS["addCellInfoCellType"]> {
  return {
    payload,
    type: "addCellInfoCellType",
  };
}

export function loadStateFromURL(
  payload: LoadStateFromURLPayload
): GetActionTypeOfReducer<typeof REDUCERS["loadStateFromURL"]> {
  return {
    payload,
    type: "loadStateFromURL",
  };
}

type GetActionTypeOfReducer<T> = T extends (
  state: never,
  action: infer Action
) => unknown
  ? Action
  : never;
