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
} from "src/views/WheresMyGene/common/store";
import {
  addGeneInfoGene,
  clearGeneInfoGene,
  closeRightSidebar,
} from "src/views/WheresMyGene/common/store/actions";
import {
  GeneExpressionSummary,
  ChartProps,
} from "src/views/WheresMyGene/common/types";

// eslint-disable-next-line sonarjs/cognitive-complexity
export const useConnect = () => {
  /* **************************************
   * Global state of the app
   * **************************************/
  const state = useContext(StateContext);
  const dispatch = useContext(DispatchContext);
  const {
    selectedGenes,
    sortBy,
    geneInfoGene,
    cellInfoCellType,
    filteredCellTypes,
    expandedTissueIds,
  } = state;
  const selectedOrganismId = state.selectedOrganismId || "";

  /* **************************************
   * Data fetching
   * **************************************/
  /**
   * Primary filter dimensions
   */
  const { data: { tissues } = {} } = usePrimaryFilterDimensions(2); //temp explicit version 2
  /**
   * Tissues keyed by ID
   */
  let tissuesByID;
  if (tissues) {
    tissuesByID = generateTermsByKey(tissues, "id");
  }
  /**
   * CellTypes keyed by Tissue Name.
   */
  const {
    data: rawCellTypesByTissueName,
    isLoading: isLoadingCellTypesByTissueName,
  } = useCellTypesByTissueName(2);
  /**
   * Gene expression summaries keyed by tissue name.
   */
  const { data: rawGeneExpressionSummariesByTissueName, isLoading } =
    useGeneExpressionSummariesByTissueName(2);

  /* ***************************************
   * Local state of the app
   * *************************************/
  /**
   * Gene expression summaries keyed by tissue name.
   */
  const [
    geneExpressionSummariesByTissueName,
    setGeneExpressionSummariesByTissueName,
  ] = useState<GeneExpressionSummariesByTissueName>(EMPTY_OBJECT);
  /**
   * Sidebar width
   */
  const [sidebarWidth, setSidebarWidth] = useState(
    FILTERS_PANEL_EXPANDED_WIDTH_PX
  );
  /**
   * Chart props for all tissues
   */
  const [allChartProps, setAllChartProps] = useState<{
    [tissue: string]: ChartProps;
  }>({});
  /**
   * Available filters - set inside Filter Panel
   */
  const [availableFilters, setAvailableFilters] =
    useState<Partial<FilterDimensions>>(EMPTY_OBJECT);
  /**
   * Whether the heatmap is scaled or not
   */
  const [isScaled, setIsScaled] = useState(true);
  /**
   * This is set in HeatMap and the value is used as
   * a list of tissues in SaveExport
   */
  const [tissuesByName, setTissuesByName] = useState<{
    [name: string]: OntologyTerm;
  }>({});
  /**
   * Cell types keyed by tissue name
   */
  const [cellTypesByTissueName, setCellTypesByTissueName] =
    useState<CellTypeByTissueName>(EMPTY_OBJECT);
  /**
   * Whether the source dataset sidebar is open or not
   */
  const [isSourceDatasetSidebarOpen, setSourceDatasetSidebarOpen] =
    useState(false);
  /**
   * Download status for SaveExport
   */
  const [downloadStatus, setDownloadStatus] = useState<{
    isLoading: boolean;
  }>({
    isLoading: false,
  });
  /**
   * Echarts renderer mode (canvas or svg) used in SaveExport
   */
  const [echartsRendererMode, setEchartsRendererMode] = useState<
    "canvas" | "svg"
  >("canvas");

  /* ***************************************
   * Derived state of the app
   * *************************************/
  /**
   * This is needed to prevent overwriting the cellTypesByTissueName state with empty
   * object when the cell types are still loading.
   */
  useEffect(() => {
    if (isLoadingCellTypesByTissueName) return;
    setCellTypesByTissueName(rawCellTypesByTissueName);
  }, [rawCellTypesByTissueName, isLoadingCellTypesByTissueName]);

  useEffect(() => {
    if (isLoading) return;

    setGeneExpressionSummariesByTissueName(
      rawGeneExpressionSummariesByTissueName
    );
  }, [rawGeneExpressionSummariesByTissueName, isLoading]);

  // TODO(thuang): Fix this complexity
  // eslint-disable-next-line sonarjs/cognitive-complexity
  const { scaledMeanExpressionMax, scaledMeanExpressionMin } = useMemo(() => {
    let min = Infinity;
    let max = -Infinity;

    for (const [tissueName, tissueSelectedCellTypes] of Object.entries(
      cellTypesByTissueName
    )) {
      const tissueSelectedCellTypeIds = tissueSelectedCellTypes.map(
        (cellType) => cellType.viewId
      );
      const tissueGeneExpressionSummaries =
        geneExpressionSummariesByTissueName[tissueName];

      if (!tissueGeneExpressionSummaries) {
        continue;
      }

      for (const selectedGeneName of selectedGenes) {
        const geneExpressionSummary =
          tissueGeneExpressionSummaries[selectedGeneName];

        if (geneExpressionSummary) {
          const { cellTypeGeneExpressionSummaries } = geneExpressionSummary;

          for (const cellTypeGeneExpressionSummary of cellTypeGeneExpressionSummaries) {
            if (
              !tissueSelectedCellTypeIds.includes(
                cellTypeGeneExpressionSummary.viewId
              )
            ) {
              continue;
            }

            const { meanExpression } = cellTypeGeneExpressionSummary;

            min = Math.min(min, meanExpression);
            max = Math.max(max, meanExpression);
          }
        }
      }
    }

    return {
      scaledMeanExpressionMax: max,
      scaledMeanExpressionMin: min,
    };
  }, [
    geneExpressionSummariesByTissueName,
    cellTypesByTissueName,
    selectedGenes,
  ]);

  /**
   * Returns an object containing gene expression summaries for each tissue name
   * based on the selected genes and gene expression summaries by tissue name.
   * If there is no expression data for a selected gene in a tissue, an empty
   * object is returned to ensure the heatmap's gene column is available.
   */
  const selectedGeneExpressionSummariesByTissueName = useMemo(() => {
    const result: { [tissueName: string]: GeneExpressionSummary[] } = {};
    for (const tissueName of Object.keys(cellTypesByTissueName)) {
      const tissueGeneExpressionSummaries =
        geneExpressionSummariesByTissueName[tissueName];

      if (!tissueGeneExpressionSummaries) continue;

      result[tissueName] = selectedGenes.map((geneName) => {
        if (tissueGeneExpressionSummaries[geneName])
          return tissueGeneExpressionSummaries[geneName];
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

  /* ***************************************
   * Handlers
   * ***************************************/
  const handleCloseRightSideBar = () => {
    if (!dispatch) return;
    dispatch(closeRightSidebar());
  };

  const handleCloseGeneInfoSideBar = () => {
    if (!dispatch) return;
    dispatch(clearGeneInfoGene());
  };

  const handleSourceDatasetButtonClick = useCallback(() => {
    setSourceDatasetSidebarOpen(!isSourceDatasetSidebarOpen);
  }, [isSourceDatasetSidebarOpen]);

  const generateGeneInfo = (gene: string) => {
    if (!dispatch) return;
    dispatch(addGeneInfoGene(gene));
  };

  return {
    sortBy,
    geneInfoGene,
    cellInfoCellType,
    filteredCellTypes,
    echartsRendererMode,
    expandedTissueIds,
    tissuesByID,
    availableFilters,
    sidebarWidth,
    cellTypesByTissueName,
    downloadStatus,
    allChartProps,
    tissuesByName,
    scaledMeanExpressionMax,
    scaledMeanExpressionMin,
    check: {
      loading: isLoading,
      scaled: isScaled,
      sourceDatasetSidebarOpen: isSourceDatasetSidebarOpen,
    },
    selected: {
      genes: selectedGenes,
      organismId: selectedOrganismId,
      geneExpressionSummariesByTissueName:
        selectedGeneExpressionSummariesByTissueName,
    },
    handle: {
      closeRightSideBar: handleCloseRightSideBar,
      closeGeneInfoSideBar: handleCloseGeneInfoSideBar,
      sourceDatasetButtonClick: handleSourceDatasetButtonClick,
      generateGeneInfo,
    },
    set: {
      sidebarWidth: setSidebarWidth,
      availableFilters: setAvailableFilters,
      tissuesByName: setTissuesByName,
      isScaled: setIsScaled,
      downloadStatus: setDownloadStatus,
      echartsRendererMode: setEchartsRendererMode,
      allChartProps: setAllChartProps,
    },
  };
};
