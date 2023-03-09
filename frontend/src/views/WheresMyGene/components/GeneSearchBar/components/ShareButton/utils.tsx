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

export const generateAndCopyShareUrl = ({
  filters,
  tissues,
  genes,
  compare,
}: {
  filters: State["selectedFilters"];
  tissues: State["selectedTissues"];
  genes: State["selectedGenes"];
  compare: State["compare"];
}) => {
  // Create a URL that contains the selected filters, tissues, and genes as params in the URL
  // This URL can be shared with others to reproduce the same view
  const url = new URL(window.location.href);
  Object.entries(stripEmptyFilters(filters)).forEach(([key, value]) => {
    url.searchParams.set(key, value.join(","));
  });
  url.searchParams.set("tissues", tissues.join(","));
  url.searchParams.set("genes", genes.join(","));
  url.searchParams.set("ver", "2");

  if (compare) {
    url.searchParams.set("compare", compare);
  }

  // Copy the URL to the clipboard
  navigator.clipboard.writeText(String(url));
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

  // Check for version
  let version = params.get("ver");
  if (!version) version = "1";
  const delimiter = version === "1" ? "-" : ",";
  // Check for tissues
  const newSelectedTissues = params.get("tissues")?.split(delimiter) || [];
  if (newSelectedTissues.length > 0) paramsToRemove.push("tissues");

  //Check for genes
  const newSelectedGenes = params.get("genes")?.split(delimiter) || [];
  if (newSelectedGenes.length > 0) paramsToRemove.push("genes");
  if (params.get("ver")) paramsToRemove.push("ver");
  if (
    Object.values(Object.keys(newSelectedFilters)).length === 0 &&
    newSelectedTissues.length === 0 &&
    newSelectedGenes.length === 0
  ) {
    return null;
  }

  // check for compare
  const newCompare = (params.get("compare") as CompareId) || undefined;

  if (newCompare) {
    paramsToRemove.push("compare");
  }

  removeParams(paramsToRemove);

  dispatch(
    loadStateFromURL({
      compare: newCompare,
      filters: newSelectedFilters,
      genes: newSelectedGenes,
      tissues: newSelectedTissues,
    })
  );

  return {
    compare: newCompare,
    filters: newSelectedFilters,
    genes: newSelectedGenes,
    tissues: newSelectedTissues,
  };
};
