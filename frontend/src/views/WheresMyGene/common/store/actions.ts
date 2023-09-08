/**
 * Create a corresponding action creator function for each action type defined
 * in `./reducer.ts`.
 */

import {
  AddCellInfoCellTypePayload,
  LoadStateFromURLPayload,
  REDUCERS,
  State,
} from "./reducer";

export function deleteSingleGene(
  geneToDelete: string
): GetActionTypeOfReducer<(typeof REDUCERS)["deleteSingleGene"]> {
  return {
    payload: geneToDelete,
    type: "deleteSingleGene",
  };
}

export function deleteSelectedGenes(): GetActionTypeOfReducer<
  (typeof REDUCERS)["deleteSelectedGenes"]
> {
  return {
    payload: null,
    type: "deleteSelectedGenes",
  };
}

export function deleteAllGenes(): GetActionTypeOfReducer<
  (typeof REDUCERS)["deleteAllGenes"]
> {
  return {
    payload: null,
    type: "deleteAllGenes",
  };
}

export function toggleGeneToDelete(
  geneToDelete: string
): GetActionTypeOfReducer<(typeof REDUCERS)["toggleGeneToDelete"]> {
  return {
    payload: geneToDelete,
    type: "toggleGeneToDelete",
  };
}

export function selectOrganism(
  organismId: State["selectedOrganismId"]
): GetActionTypeOfReducer<(typeof REDUCERS)["selectOrganism"]> {
  return {
    payload: organismId,
    type: "selectOrganism",
  };
}

export function selectGenes(
  genes: State["selectedGenes"]
): GetActionTypeOfReducer<(typeof REDUCERS)["selectGenes"]> {
  return {
    payload: genes,
    type: "selectGenes",
  };
}

export function addSelectedGenes(
  genes: State["selectedGenes"]
): GetActionTypeOfReducer<(typeof REDUCERS)["addSelectedGenes"]> {
  return {
    payload: genes,
    type: "addSelectedGenes",
  };
}

export function selectTissues(
  tissues: State["selectedTissues"]
): GetActionTypeOfReducer<(typeof REDUCERS)["selectTissues"]> {
  return {
    payload: tissues,
    type: "selectTissues",
  };
}

export function selectSortBy(
  sortBy: Partial<State["sortBy"]>
): GetActionTypeOfReducer<(typeof REDUCERS)["selectSortBy"]> {
  return {
    payload: sortBy,
    type: "selectSortBy",
  };
}

export function selectCompare(
  Compare: State["compare"]
): GetActionTypeOfReducer<(typeof REDUCERS)["selectCompare"]> {
  return {
    payload: Compare,
    type: "selectCompare",
  };
}

export function resetGenesToDelete(): GetActionTypeOfReducer<
  (typeof REDUCERS)["resetGenesToDelete"]
> {
  return {
    payload: null,
    type: "resetGenesToDelete",
  };
}

export function selectFilters(
  key: keyof State["selectedFilters"],
  options: string[]
): GetActionTypeOfReducer<(typeof REDUCERS)["selectFilters"]> {
  return {
    payload: { key, options },
    type: "selectFilters",
  };
}

export function setSnapshotId(
  snapshotId: State["snapshotId"]
): GetActionTypeOfReducer<(typeof REDUCERS)["setSnapshotId"]> {
  return {
    payload: snapshotId,
    type: "setSnapshotId",
  };
}

export function addCellInfoCellType(
  payload: AddCellInfoCellTypePayload
): GetActionTypeOfReducer<(typeof REDUCERS)["addCellInfoCellType"]> {
  return {
    payload,
    type: "addCellInfoCellType",
  };
}

export function addGeneInfoGene(
  payload: string
): GetActionTypeOfReducer<(typeof REDUCERS)["addGeneInfoGene"]> {
  return {
    payload,
    type: "addGeneInfoGene",
  };
}

export function clearGeneInfoGene(): GetActionTypeOfReducer<
  (typeof REDUCERS)["clearGeneInfoGene"]
> {
  return {
    payload: null,
    type: "clearGeneInfoGene",
  };
}

export function clearCellInfoCellType(): GetActionTypeOfReducer<
  (typeof REDUCERS)["clearCellInfoCellType"]
> {
  return {
    payload: null,
    type: "clearCellInfoCellType",
  };
}

export function closeRightSidebar(): GetActionTypeOfReducer<
  (typeof REDUCERS)["closeRightSidebar"]
> {
  return {
    payload: null,
    type: "closeRightSidebar",
  };
}

export function selectGeneInfoFromXAxis(
  payload: string
): GetActionTypeOfReducer<(typeof REDUCERS)["selectGeneInfoFromXAxis"]> {
  return {
    payload,
    type: "selectGeneInfoFromXAxis",
  };
}

export function loadStateFromURL(
  payload: LoadStateFromURLPayload
): GetActionTypeOfReducer<(typeof REDUCERS)["loadStateFromURL"]> {
  return {
    payload,
    type: "loadStateFromURL",
  };
}

export function setXAxisHeight(
  payload: number
): GetActionTypeOfReducer<(typeof REDUCERS)["setXAxisHeight"]> {
  return {
    payload,
    type: "setXAxisHeight",
  };
}

interface SetFilteredCellTypesPayload {
  filteredCellTypes: State["filteredCellTypes"];
  filteredCellTypeIds: State["filteredCellTypeIds"];
}

export function setFilteredCellTypes(
  payload: SetFilteredCellTypesPayload
): GetActionTypeOfReducer<(typeof REDUCERS)["setFilteredCellTypes"]> {
  return {
    payload,
    type: "setFilteredCellTypes",
  };
}

type GetActionTypeOfReducer<T> = T extends (
  state: never,
  action: infer Action
) => unknown
  ? Action
  : never;
