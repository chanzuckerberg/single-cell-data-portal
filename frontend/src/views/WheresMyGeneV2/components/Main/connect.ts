import { useCallback, useContext, useEffect, useMemo, useState } from "react";

import { EMPTY_ARRAY, EMPTY_OBJECT } from "src/common/constants/utils";

import {
  CellTypeByTissueName,
  FilterDimensions,
  GeneExpressionSummariesByTissueName,
  generateTermsByKey,
  OntologyTerm,
  useCellTypesByTissueName,
  useGeneExpressionSummariesByTissueName,
  usePrimaryFilterDimensions,
} from "src/common/queries/wheresMyGene";

import { FILTERS_PANEL_EXPANDED_WIDTH_PX } from "src/components/common/SideBar";

import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGeneV2/common/store";

import {
  GeneExpressionSummary,
  ChartProps,
} from "src/views/WheresMyGeneV2/common/types";
import { useScaledMeanExpression } from "src/views/WheresMyGeneV2/components/Main/hooks";

export function useConnect() {
  const state = useContext(StateContext);
  const dispatch = useContext(DispatchContext);

  const {
    selectedGenes,
    sortBy,
    geneInfoGene,
    cellInfoCellType,
    filteredCellTypes,
    expandedTissueIds,
    selectedFilters: { tissues: filteredTissueIds },
  } = state;

  const selectedOrganismId = state.selectedOrganismId || "";

  const { data: { tissues } = {} } = usePrimaryFilterDimensions();

  let tissuesByID;

  if (tissues) {
    tissuesByID = generateTermsByKey(tissues, "id");
  }

  const [allChartProps, setAllChartProps] = useState<{
    [tissue: string]: ChartProps;
  }>({});

  const [availableFilters, setAvailableFilters] =
    useState<Partial<FilterDimensions>>(EMPTY_OBJECT);

  const [isScaled, setIsScaled] = useState(true);

  // This is set in HeatMap and the value is used as a list of tissues in SaveExport
  const [tissuesByName, setTissuesByName] = useState<{
    [name: string]: OntologyTerm;
  }>({});

  //(seve): These useEffects are deceptively simple.
  // Their purpose is to avoid updating the state with null/empty values while we're waiting for the api to return data.

  const {
    data: rawCellTypesByTissueName,
    isLoading: isLoadingCellTypesByTissueName,
  } = useCellTypesByTissueName();

  const [cellTypesByTissueName, setCellTypesByTissueName] =
    useState<CellTypeByTissueName>(EMPTY_OBJECT);

  // This is needed to prevent overwriting the cellTypesByTissueName state with empty
  useEffect(() => {
    if (isLoadingCellTypesByTissueName) return;

    setCellTypesByTissueName(rawCellTypesByTissueName);
  }, [rawCellTypesByTissueName, isLoadingCellTypesByTissueName]);

  /**
   * This holds ALL the geneData we have loaded from the API, including previously
   * and currently selected genes.
   * We use `selectedGeneData` to subset the data to only the genes that are
   * currently selected.
   */
  const { data: rawGeneExpressionSummariesByTissueName, isLoading } =
    useGeneExpressionSummariesByTissueName();

  const [
    geneExpressionSummariesByTissueName,
    setGeneExpressionSummariesByTissueName,
  ] = useState<GeneExpressionSummariesByTissueName>(EMPTY_OBJECT);

  useEffect(() => {
    if (isLoading) return;

    setGeneExpressionSummariesByTissueName(
      rawGeneExpressionSummariesByTissueName
    );
  }, [rawGeneExpressionSummariesByTissueName, isLoading]);

  const { scaledMeanExpressionMax, scaledMeanExpressionMin } =
    useScaledMeanExpression({
      cellTypesByTissueName,
      geneExpressionSummariesByTissueName,
      selectedGenes,
      expandedTissueIds,
      filteredCellTypes,
      filteredTissueIds,
    });

  const selectedGeneExpressionSummariesByTissueName = useMemo(() => {
    const result: { [tissueName: string]: GeneExpressionSummary[] } = {};

    for (const tissueName of Object.keys(cellTypesByTissueName)) {
      const tissueGeneExpressionSummaries =
        geneExpressionSummariesByTissueName[tissueName];

      if (!tissueGeneExpressionSummaries) continue;

      result[tissueName] = selectedGenes.map((geneName) => {
        // early return to avoid unnecessary generation of empty object
        if (tissueGeneExpressionSummaries[geneName])
          return tissueGeneExpressionSummaries[geneName];

        // (thuang): This is needed to ensure the heatmap's gene column
        // is available even if there's no expression data for the column.
        // Otherwise the heatmap columns and column labels won't match up
        // where there's holes in the data.
        return {
          cellTypeGeneExpressionSummaries: EMPTY_ARRAY,
          name: geneName,
        };
      });
    }

    return result;
  }, [
    geneExpressionSummariesByTissueName,
    selectedGenes,
    cellTypesByTissueName,
  ]);
  const [isSourceDatasetSidebarOpen, setSourceDatasetSidebarOpen] =
    useState(false);

  const handleSourceDatasetButtonClick = useCallback(() => {
    setSourceDatasetSidebarOpen(!isSourceDatasetSidebarOpen);
  }, [isSourceDatasetSidebarOpen]);

  const [downloadStatus, setDownloadStatus] = useState<{
    isLoading: boolean;
  }>({
    isLoading: false,
  });

  const [echartsRendererMode, setEchartsRendererMode] = useState<
    "canvas" | "svg"
  >("canvas");

  // Sidebar width
  const [sidebarWidth, setSidebarWidth] = useState(
    FILTERS_PANEL_EXPANDED_WIDTH_PX
  );

  return {
    dispatch,
    sidebarWidth,
    setSidebarWidth,
    isLoading,
    availableFilters,
    setAvailableFilters,
    isScaled,
    setIsScaled,
    handleSourceDatasetButtonClick,
    downloadStatus,
    setDownloadStatus,
    echartsRendererMode,
    setEchartsRendererMode,
    allChartProps,
    setAllChartProps,
    tissuesByName,
    setTissuesByName,
    sortBy,
    geneInfoGene,
    cellInfoCellType,
    tissuesByID,
    cellTypesByTissueName,
    selectedGenes,
    expandedTissueIds,
    filteredCellTypes,
    scaledMeanExpressionMax,
    scaledMeanExpressionMin,
    isSourceDatasetSidebarOpen,
    selectedOrganismId,
    selectedGeneExpressionSummariesByTissueName,
  };
}
