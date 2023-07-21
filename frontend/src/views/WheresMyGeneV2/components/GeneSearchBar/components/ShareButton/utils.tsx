import { Dispatch } from "react";
import { isSSR } from "src/common/utils/isSSR";
import { removeParams } from "src/common/utils/removeParams";
import { CompareId } from "src/views/WheresMyGene/common/constants";
import { State } from "src/views/WheresMyGene/common/store";
import { loadStateFromURL } from "src/views/WheresMyGene/common/store/actions";
import {
  LoadStateFromURLPayload,
  PayloadAction,
} from "src/views/WheresMyGene/common/store/reducer";

const HUMAN_ORGANISM_ID = "NCBITaxon:9606";

export const LATEST_SHARE_LINK_VERSION = "2";

export const generateAndCopyShareUrl = ({
  filters,
  organism,
  genes,
  compare,
  cellTypes = [],
  copyToClipboard = true,
}: {
  filters: State["selectedFilters"];
  organism: State["selectedOrganismId"];
  genes: State["selectedGenes"];
  compare: State["compare"];
  cellTypes?: State["filteredCellTypes"];
  copyToClipboard?: boolean;
}) => {
  // Create a URL that contains the selected filters, tissues, and genes as params in the URL
  // This URL can be shared with others to reproduce the same view
  const url = new URL(window.location.href);
  // human is empty default
  if (organism && organism !== HUMAN_ORGANISM_ID) {
    url.searchParams.set("organism", organism);
  }
  Object.entries(stripEmptyFilters(filters)).forEach(([key, value]) => {
    url.searchParams.set(key, value.join(","));
  });
  url.searchParams.set("genes", genes.join(","));
  url.searchParams.set("ver", LATEST_SHARE_LINK_VERSION);

  if (compare) {
    url.searchParams.set("compare", compare);
  }

  if (cellTypes.length > 0) {
    url.searchParams.set("cellTypes", cellTypes.join(","));
  }

  const urlString = String(url);

  if (copyToClipboard) {
    // Copy the URL to the clipboard
    navigator.clipboard.writeText(urlString);
  }

  return urlString;
};

const stripEmptyFilters = (
  filters: State["selectedFilters"]
): Partial<State["selectedFilters"]> => {
  // Remove filters that don't have any values selected
  const strippedFilters: Partial<State["selectedFilters"]> = {};
  Object.entries(filters).forEach(([key, value]) => {
    if (value.length > 0) {
      strippedFilters[key as keyof State["selectedFilters"]] = value;
    }
  });
  return strippedFilters;
};

export const loadStateFromQueryParams = (
  params: URLSearchParams,
  selectedFilters: State["selectedFilters"],
  dispatch: Dispatch<PayloadAction<LoadStateFromURLPayload>>
): LoadStateFromURLPayload | null => {
  if (isSSR()) return null;

  const paramsToRemove = [];

  // Checks if the URL has any query params related to shared filter state
  // If so, it will update the state with the values from the URL

  // Check for version
  const version = params.get("ver") || "1";
  if (params.get("ver")) paramsToRemove.push("ver");

  // delimiter changed from - to , in version 2
  const delimiter = version > "1" ? "," : "-";

  // Check for filter properties

  const newSelectedFilters: Partial<State["selectedFilters"]> = {};
  Object.keys(selectedFilters).forEach((key) => {
    const value = params.get(key);
    if (value) {
      newSelectedFilters[key as keyof State["selectedFilters"]] =
        value.split(delimiter);
      paramsToRemove.push(key);
    }
  });

  //Check for organism
  const newSelectedOrganism = params.get("organism") || HUMAN_ORGANISM_ID;
  if (newSelectedOrganism) {
    paramsToRemove.push("organism");
  }

  //Check for genes
  const newSelectedGenes = params.get("genes")?.split(delimiter) || [];
  if (newSelectedGenes.length > 0) paramsToRemove.push("genes");

  //Check for cell types
  const newFilteredCellTypes = params.get("cellTypes")?.split(delimiter) || [];
  if (newFilteredCellTypes.length > 0) paramsToRemove.push("cellTypes");

  removeParams(paramsToRemove);

  // If there are no filters and genes selected, don't update the state
  if (
    Object.values(Object.keys(newSelectedFilters)).length === 0 &&
    newSelectedGenes.length === 0
  ) {
    return null;
  }

  // check for compare
  const newCompare = (params.get("compare") as CompareId) || undefined;

  if (newCompare) {
    paramsToRemove.push("compare");
  }

  dispatch(
    loadStateFromURL({
      compare: newCompare,
      filters: newSelectedFilters,
      organism: newSelectedOrganism,
      genes: newSelectedGenes,
      cellTypes: newFilteredCellTypes,
    })
  );

  return {
    compare: newCompare,
    filters: newSelectedFilters,
    organism: newSelectedOrganism,
    genes: newSelectedGenes,
    cellTypes: newFilteredCellTypes,
  };
};
