import { Dispatch } from "react";
import { isSSR } from "src/common/utils/isSSR";
import { removeParams } from "src/common/utils/removeParams";
import { State } from "src/views/WheresMyGene/common/store";
import { loadStateFromURL } from "src/views/WheresMyGene/common/store/actions";
import {
  LoadStateFromURLPayload,
  PayloadAction,
} from "src/views/WheresMyGene/common/store/reducer";

export const generateAndCopyShareUrl = (
  filters: State["selectedFilters"],
  tissues: State["selectedTissues"],
  genes: State["selectedGenes"]
) => {
  // Create a URL that contains the selected filters, tissues, and genes as params in the URL
  // This URL can be shared with others to reproduce the same view
  const url = new URL(window.location.href);
  Object.entries(stripEmptyFilters(filters)).forEach(([key, value]) => {
    url.searchParams.set(key, value.join(","));
  });
  url.searchParams.set("tissues", tissues.join(","));
  url.searchParams.set("genes", genes.join(","));

  // Copy the URL to the clipboard
  navigator.clipboard.writeText(url.toString());
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

  // Check for filter properties

  const newSelectedFilters: Partial<State["selectedFilters"]> = {};
  Object.keys(selectedFilters).forEach((key) => {
    const value = params.get(key);
    if (value) {
      newSelectedFilters[key as keyof State["selectedFilters"]] =
        value.split("-");
      paramsToRemove.push(key);
    }
  });

  // Check for tissues
  const newSelectedTissues = params.get("tissues")?.split(",") || [];
  if (newSelectedTissues.length > 0) paramsToRemove.push("tissues");

  //Check for genes
  const newSelectedGenes = params.get("genes")?.split(",") || [];
  if (newSelectedGenes.length > 0) paramsToRemove.push("genes");

  if (
    Object.values(Object.keys(newSelectedFilters)).length === 0 &&
    newSelectedTissues.length === 0 &&
    newSelectedGenes.length === 0
  )
    return null;

  removeParams(paramsToRemove);

  dispatch(
    loadStateFromURL({
      filters: newSelectedFilters,
      tissues: newSelectedTissues,
      genes: newSelectedGenes,
    })
  );
  return {
    filters: newSelectedFilters,
    tissues: newSelectedTissues,
    genes: newSelectedGenes,
  };
};
