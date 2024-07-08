import { NextRouter } from "next/router";
import { Dispatch } from "react";
import { ROUTES } from "src/common/constants/routes";
import { TissueMetadataQueryResponse } from "src/common/queries/cellGuide";
import { isSSR } from "src/common/utils/isSSR";
import { removeParams } from "src/common/utils/removeParams";
import { CompareId } from "src/views/WheresMyGeneV2/common/constants";
import { State } from "src/views/WheresMyGeneV2/common/store";
import { loadStateFromURL } from "src/views/WheresMyGeneV2/common/store/actions";
import {
  LoadStateFromURLPayload,
  PayloadAction,
} from "src/views/WheresMyGeneV2/common/store/reducer";
import { CellType } from "src/views/WheresMyGeneV2/common/types";

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

  if (compare) {
    url.searchParams.set("compare", compare);
  }

  Object.entries(stripEmptyFilters(filters)).forEach(([key, value]) => {
    url.searchParams.set(key, value.join(","));
  });
  url.searchParams.set("genes", genes.join(","));

  if (cellTypes.length > 0) {
    url.searchParams.set("cellTypes", cellTypes.join(","));
  }

  url.searchParams.set("ver", LATEST_SHARE_LINK_VERSION);

  const urlString = String(url);

  if (copyToClipboard) {
    // Copy the URL to the clipboard
    navigator.clipboard.writeText(urlString);
  }

  return urlString;
};

export const generateDifferentialExpressionUrl = ({
  filters,
  organism,
  tissue,
  cellType,
}: {
  filters: State["selectedFilters"];
  organism: State["selectedOrganismId"];
  tissue: string;
  cellType: string;
}) => {
  // Create a URL that contains the selected filters, cell type, and tissue as params in the URL
  // This URL can be shared with others to reproduce the same view
  const url = new URL(ROUTES.DE, window.location.origin);

  // human is empty default
  if (organism && organism !== HUMAN_ORGANISM_ID) {
    url.searchParams.set("organism", organism);
  }

  Object.entries(stripEmptyFilters(filters)).forEach(([key, value]) => {
    if (["diseases", "ethnicities", "publications", "sexes"].includes(key)) {
      url.searchParams.set(key, value.join(","));
    }
  });

  url.searchParams.set("celltypes", cellType);
  url.searchParams.set("tissues", tissue);
  return String(url);
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

export const loadStateFromQueryParams = ({
  cellTypesByName,
  dispatch,
  params,
  router,
  selectedFilters,
  tissues,
}: {
  params: URLSearchParams;
  selectedFilters: State["selectedFilters"];
  dispatch: Dispatch<PayloadAction<LoadStateFromURLPayload>>;
  tissues: TissueMetadataQueryResponse;
  cellTypesByName: { [name: string]: CellType };
  router: NextRouter;
}): LoadStateFromURLPayload | null => {
  if (isSSR()) return null;

  const paramsToRemove: string[] = [];

  // Checks if the URL has any query params related to shared filter state
  // If so, it will update the state with the values from the URL

  const { delimiter } = getMetadataForShareLinkParsing({
    params,
    paramsToRemove,
  });

  // Check for filter properties
  const newSelectedFilters = getNewSelectedFilters({
    tissues,
    params,
    paramsToRemove,
    selectedFilters,
    delimiter,
  });

  // Check for genes
  const newSelectedGenes = params.get("genes")?.split(delimiter) || [];
  if (newSelectedGenes.length > 0) paramsToRemove.push("genes");

  // If there are no filters and genes selected, don't update the state
  if (
    Object.values(Object.keys(newSelectedFilters)).length === 0 &&
    newSelectedGenes.length === 0
  ) {
    return null;
  }

  // Check for organism
  const newSelectedOrganism = params.get("organism") || HUMAN_ORGANISM_ID;

  if (newSelectedOrganism) {
    paramsToRemove.push("organism");
  }

  const cellTypeNames = Object.keys(cellTypesByName);

  // Check for cell types
  const newFilteredCellTypes =
    params
      .get("cellTypes")
      ?.split(delimiter)
      .filter((cellType) => {
        if (cellTypeNames.indexOf(cellType) !== -1) return true;
        // Pop toast here
        console.warn(`Cell type ${cellType} not found in cell types`);
        return false;
      }) || [];

  const newFilteredCellTypeIds = newFilteredCellTypes.map(
    (cellTypeName) => cellTypesByName[cellTypeName].id
  );

  if (newFilteredCellTypes.length > 0) paramsToRemove.push("cellTypes");

  // Check for compare
  const newCompare = (params.get("compare") as CompareId) || undefined;

  if (newCompare) {
    paramsToRemove.push("compare");
  }

  /**
   * (thuang): Please makes sure we only remove params AFTER pushing all params
   * to `paramsToRemove`
   */
  removeParams({ params: paramsToRemove, router });

  const payload = {
    compare: newCompare,
    filters: newSelectedFilters,
    organism: newSelectedOrganism,
    genes: newSelectedGenes,
    cellTypes: newFilteredCellTypes,
    cellTypeIds: newFilteredCellTypeIds,
  };

  dispatch(loadStateFromURL(payload));

  return payload;
};

function getMetadataForShareLinkParsing({
  params,
  paramsToRemove,
}: {
  params: URLSearchParams;
  paramsToRemove: string[];
}) {
  const version = params.get("ver") || "1";
  if (params.get("ver")) paramsToRemove.push("ver");

  return {
    // delimiter changed from `-` to `,` in version 2
    delimiter: version > "1" ? "," : "-",
  };
}

function getNewSelectedFilters({
  tissues,
  params,
  paramsToRemove,
  selectedFilters,
  delimiter,
}: {
  tissues: TissueMetadataQueryResponse;
  params: URLSearchParams;
  paramsToRemove: string[];
  selectedFilters: State["selectedFilters"];
  delimiter: string;
}) {
  // Check for filter properties
  const newSelectedFilters: Partial<State["selectedFilters"]> = {};

  /**
   * (thuang): Map tissue names to IDs for backwards compatibility
   */
  const tissueIdsByName = new Map(
    Object.entries(tissues ?? {}).map(([id, tissue]) => [tissue.name, id])
  );

  /**
   * (thuang): Warning - Map doesn't work with Object.keys(), so we need to use
   * Map.keys() instead
   */
  const allTissueNames = Array.from(tissueIdsByName.keys());
  const allTissueIds = Object.keys(tissues ?? {});

  Object.keys(selectedFilters).forEach((key) => {
    const value = params.get(key);

    if (!value) return;

    /**
     * (thuang): Special case to allow tissues to be specified by name or ID
     * https://github.com/chanzuckerberg/single-cell-data-portal/issues/5358
     */
    if (key === "tissues") {
      // Check if the values are tissue names or IDs
      const tissueParams = value.split(delimiter);
      const tissueIds = [];
      const tissueNames = [];

      for (const tissueParam of tissueParams) {
        if (
          tissueParam.includes("UBERON:") &&
          allTissueIds.indexOf(tissueParam) !== -1
        ) {
          tissueIds.push(tissueParam);
        } else if (allTissueNames.indexOf(tissueParam) !== -1) {
          tissueNames.push(tissueParam);
        }
      }

      newSelectedFilters["tissues"] = [
        ...tissueIds,
        ...(tissueNames.map((tissueName) =>
          tissueIdsByName.get(tissueName)
        ) as string[]),
      ];
    } else {
      // TODO(seve): add input validation for all filter types, ensure that the values are present in the data
      newSelectedFilters[key as keyof State["selectedFilters"]] =
        value.split(delimiter);
    }

    paramsToRemove.push(key);
  });

  return newSelectedFilters;
}
