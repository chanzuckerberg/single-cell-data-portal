import { DrawerSize } from "@blueprintjs/core";
import Head from "next/head";
import { useCallback, useContext, useEffect, useMemo, useState } from "react";
import { EMPTY_ARRAY, EMPTY_OBJECT } from "src/common/constants/utils";
import {
  CellTypeByTissueName,
  FilterDimensions,
  GeneExpressionSummariesByTissueName,
  generateTermsByKey,
  getOntologyTermIdFromCellTypeViewId,
  OntologyTerm,
  useCellTypesByTissueName,
  useGeneExpressionSummariesByTissueName,
  usePrimaryFilterDimensions,
} from "src/common/queries/wheresMyGene";
import SideBar, {
  FILTERS_PANEL_EXPANDED_WIDTH_PX,
} from "src/components/common/SideBar";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGeneV2/common/store";
import {
  addGeneInfoGene,
  clearGeneInfoGene,
  closeRightSidebar,
} from "src/views/WheresMyGeneV2/common/store/actions";
import {
  GeneExpressionSummary,
  ChartProps,
} from "src/views/WheresMyGeneV2/common/types";
import CellInfoSideBar from "src/views/WheresMyGeneV2/components/CellInfoSideBar";
import Filters from "src/views/WheresMyGeneV2/components/Filters";
import GeneInfoSideBar from "src/components/GeneInfoSideBar";

import InfoPanel from "src/views/WheresMyGeneV2/components/InfoPanel";
import Legend from "src/views/WheresMyGeneV2/components/InfoPanel/components/Legend";
import Loader from "src/views/WheresMyGeneV2/components/Loader";
import {
  StyledBannerContainer,
  StyledSidebarDrawer,
} from "src/views/WheresMyGene/components/Main/style";
import ScreenTint from "src/views/WheresMyGeneV2/components/ScreenTint";
import {
  SideBarPositioner,
  SideBarWrapper,
  Top,
  Wrapper,
} from "src/views/WheresMyGene/style";
import { View } from "src/views/globalStyle";
import HeatMap from "../HeatMap";
import BottomBanner from "src/components/BottomBanner";
import { CELL_INFO_SIDEBAR_WIDTH_PX } from "src/views/WheresMyGeneV2/components/CellInfoSideBar/style";
import { UnderlyingDataChangeBanner } from "../GeneSearchBar/components/SaveExport/ExportBanner";
import { GENE_EXPRESSION_BANNER_SURVEY_LINK } from "src/common/constants/airtableLinks";
import { StyledRightSideBar } from "./style";

export const INFO_PANEL_WIDTH_PX = 320;

export default function WheresMyGene(): JSX.Element {
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

  const { data: { tissues } = {} } = usePrimaryFilterDimensions(2); //temp explicit version 2

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
  } = useCellTypesByTissueName(2);

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
    useGeneExpressionSummariesByTissueName(2);

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

  // TODO(thuang): Fix this complexity
  // eslint-disable-next-line sonarjs/cognitive-complexity
  const { scaledMeanExpressionMax, scaledMeanExpressionMin } = useMemo(() => {
    let min = Infinity;
    let max = -Infinity;
    for (const [tissueName, tissueSelectedCellTypes] of Object.entries(
      cellTypesByTissueName
    )) {
      const tissueId = tissueSelectedCellTypes.at(-1)?.id ?? "";
      const tissueSelectedCellTypeIds = tissueSelectedCellTypes.map(
        (cellType) => cellType.viewId
      );
      const tissueSelectedCellTypeNames = tissueSelectedCellTypes.map(
        (cellType) => cellType.cellTypeName
      );
      // get object of cell type id to cell type name
      const cellTypeNameById: { [id: string]: string } = {};
      tissueSelectedCellTypes.forEach((cellType) => {
        cellTypeNameById[cellType.id] = cellType.cellTypeName;
      });
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
            // get term before $
            const cellTypeGeneExpressionSummaryId =
              getOntologyTermIdFromCellTypeViewId(
                cellTypeGeneExpressionSummary.viewId
              );
            const isCellType =
              !cellTypeGeneExpressionSummaryId.startsWith("UBERON");

            /*
            This is to dynamically set the min and max based on the data that is visible to users.
            Skip conditions:
            - If the cell type is not part of an expanded tissue
            - If the tissue does not contain any of the filtered cell types and at least one cell type is filtered
            - If the cell type is not included in the filtered cell types and at least one cell type is filtered
            - If the tissue id is not included in the filtered tissue ids and at least one tissue is filtered
            - If the cell type gene expression summary view id is not included in the available view ids
            The above conditions capture all possible scenarios where the data is not visible to users.
            */
            if (
              // If the cell type is not part of an expanded tissue
              (isCellType && !expandedTissueIds.includes(tissueId)) ||
              // If the tissue does not contain any of the filtered cell types and at least one cell type is filtered
              (!isCellType &&
                filteredCellTypes.length > 0 &&
                !tissueSelectedCellTypeNames.filter((value) =>
                  filteredCellTypes.includes(value)
                ).length) ||
              // If the cell type is not included in the filtered cell types and at least one cell type is filtered
              (isCellType &&
                filteredCellTypes.length > 0 &&
                !filteredCellTypes.includes(
                  cellTypeNameById?.[cellTypeGeneExpressionSummaryId]
                )) ||
              // If the tissue id is not included in the filtered tissue ids and at least one tissue is filtered
              (!isCellType &&
                filteredTissueIds.length > 0 &&
                !filteredTissueIds.includes(tissueId)) ||
              // If the cell type gene expression summary view id is not included in the available view ids
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
    filteredTissueIds,
    expandedTissueIds,
    filteredCellTypes,
  ]);

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

  const handleCloseRightSideBar = () => {
    if (!dispatch) return;
    dispatch(closeRightSidebar());
  };

  const handleCloseGeneInfoSideBar = () => {
    if (!dispatch) return;
    dispatch(clearGeneInfoGene());
  };

  const generateGeneInfo = (gene: string) => {
    if (!dispatch) return;
    dispatch(addGeneInfoGene(gene));
  };

  // Sidebar width
  const [sidebarWidth, setSidebarWidth] = useState(
    FILTERS_PANEL_EXPANDED_WIDTH_PX
  );

  return (
    <>
      <Head>
        <title>Gene Expression - CZ CELLxGENE Discover</title>
      </Head>

      <SideBar
        label="Filters"
        SideBarWrapperComponent={SideBarWrapper}
        SideBarPositionerComponent={SideBarPositioner}
        testId="filters-panel"
        wmgSideBar
        onWidthChange={setSidebarWidth}
      >
        <Filters
          isLoading={isLoading}
          availableFilters={availableFilters}
          setAvailableFilters={setAvailableFilters}
          setIsScaled={setIsScaled}
        />
      </SideBar>
      {cellInfoCellType && tissuesByID ? (
        <StyledRightSideBar width={CELL_INFO_SIDEBAR_WIDTH_PX}>
          <CellInfoSideBar
            generateGeneInfo={generateGeneInfo}
            cellInfoCellType={cellInfoCellType}
            tissueInfo={tissuesByID[cellInfoCellType.tissueID]}
            handleClose={handleCloseRightSideBar}
            title={`${cellInfoCellType.cellType.name}`}
          />

          {
            // Split right sidebar view if fmg AND gene info is populated
            geneInfoGene && (
              <GeneInfoSideBar
                geneInfoGene={geneInfoGene}
                handleClose={handleCloseGeneInfoSideBar}
                title={`${geneInfoGene}`}
              />
            )
          }
        </StyledRightSideBar>
      ) : (
        // Gene info full right sidebar length
        geneInfoGene && (
          <StyledRightSideBar>
            <GeneInfoSideBar
              geneInfoGene={geneInfoGene}
              handleClose={handleCloseGeneInfoSideBar}
              title={geneInfoGene}
            />
          </StyledRightSideBar>
        )
      )}

      <View id="view" overflow="hidden">
        <Wrapper>
          {isLoading && <Loader />}

          {/* Used for PNG and SVG exports to render message banner to render in output */}
          {downloadStatus.isLoading && (
            <StyledBannerContainer>
              <UnderlyingDataChangeBanner />
            </StyledBannerContainer>
          )}

          <Top id="top-legend">
            <Legend
              selectedCellTypes={cellTypesByTissueName}
              selectedGenes={selectedGenes}
              isScaled={isScaled}
              handleRightSidebarButtonClick={handleSourceDatasetButtonClick}
              setDownloadStatus={setDownloadStatus}
              setEchartsRendererMode={setEchartsRendererMode}
              allChartProps={allChartProps}
              availableFilters={availableFilters}
              tissues={tissuesByName}
              expandedTissueIds={expandedTissueIds}
              filteredCellTypes={filteredCellTypes}
              maxExpression={scaledMeanExpressionMax}
            />
          </Top>

          <StyledSidebarDrawer
            position="right"
            isOpen={isSourceDatasetSidebarOpen}
            title="Source Data"
            canEscapeKeyClose={true}
            canOutsideClickClose={true}
            onClose={handleSourceDatasetButtonClick}
            size={DrawerSize.SMALL}
          >
            <InfoPanel />
          </StyledSidebarDrawer>

          <ScreenTint isDownloading={downloadStatus} />

          <HeatMap
            /**
             * (thuang): Use the selected organism ID to reset the heatmap state
             */
            key={selectedOrganismId}
            echartsRendererMode={echartsRendererMode}
            cellTypeSortBy={sortBy.cellTypes}
            geneSortBy={sortBy.genes}
            isScaled={isScaled}
            isLoadingAPI={isLoading}
            cellTypes={cellTypesByTissueName}
            genes={selectedGenes}
            selectedGeneExpressionSummariesByTissueName={
              selectedGeneExpressionSummariesByTissueName
            }
            scaledMeanExpressionMax={scaledMeanExpressionMax}
            scaledMeanExpressionMin={scaledMeanExpressionMin}
            allChartProps={allChartProps}
            setAllChartProps={setAllChartProps}
            setTissuesByName={setTissuesByName}
            tissuesByName={tissuesByName}
            /**
             * (thuang): This is needed to reposition gene search bar when the
             * sidebar width expands/collapses
             */
            sidebarWidth={sidebarWidth}
          />
        </Wrapper>
      </View>

      <BottomBanner
        airtableLink={GENE_EXPRESSION_BANNER_SURVEY_LINK}
        includeSurveyLink
      />
    </>
  );
}
