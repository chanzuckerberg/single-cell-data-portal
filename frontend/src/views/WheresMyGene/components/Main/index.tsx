import { DrawerSize } from "@blueprintjs/core";
import Head from "next/head";
import React, {
  useCallback,
  useContext,
  useEffect,
  useMemo,
  useRef,
  useState,
} from "react";
import { EMPTY_ARRAY, EMPTY_OBJECT } from "src/common/constants/utils";
import {
  CellTypeByTissueName,
  FilterDimensions,
  GeneExpressionSummariesByTissueName,
  generateTermsByKey,
  useCellTypesByTissueName,
  useGeneExpressionSummariesByTissueName,
  usePrimaryFilterDimensions,
} from "src/common/queries/wheresMyGene";
import SideBar from "src/components/common/SideBar";
import { View } from "../../../globalStyle";
import { DispatchContext, StateContext } from "../../common/store";
import {
  addGeneInfoGene,
  clearGeneInfoGene,
  closeRightSidebar,
  deleteSelectedGenes,
} from "../../common/store/actions";
import { ChartProps, GeneExpressionSummary } from "../../common/types";
import { SideBarPositioner, SideBarWrapper, Top, Wrapper } from "../../style";
import CellInfoBar from "../CellInfoSideBar";
import GeneInfoBar from "../../../../components/GeneInfoSideBar";
import Filters from "../Filters";
import GeneSearchBar from "../GeneSearchBar";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../GeneSearchBar/components/SaveExport";
import GetStarted from "../GetStarted";
import HeatMap from "../HeatMap";
import InfoPanel from "../InfoPanel";
import Legend from "../InfoPanel/components/Legend";
import Loader from "../Loader";
import ScreenTint from "../ScreenTint";
import { StyledBannerContainer, StyledSidebarDrawer } from "./style";
import RightSideBar from "../../../../components/common/RightSideBar";
import { UnderlyingDataChangeBanner } from "../GeneSearchBar/components/SaveExport/ExportBanner";
import BottomBanner from "src/components/BottomBanner";
import { CELL_INFO_SIDEBAR_WIDTH_PX } from "../CellInfoSideBar/style";

export const INFO_PANEL_WIDTH_PX = 320;

export default function WheresMyGene(): JSX.Element {
  const state = useContext(StateContext);
  const dispatch = useContext(DispatchContext);

  const {
    selectedGenes,
    selectedTissues,
    sortBy,
    geneInfoGene,
    cellInfoCellType,
  } = state;

  const selectedOrganismId = state.selectedOrganismId || "";

  const { data: { tissues: allTissues } = {} } = usePrimaryFilterDimensions();

  let tissuesByID;
  if (allTissues) {
    tissuesByID = generateTermsByKey(allTissues, "id");
  }

  const [allChartProps, setAllChartProps] = useState<{
    [tissue: string]: ChartProps;
  }>({});

  const [availableFilters, setAvailableFilters] =
    useState<Partial<FilterDimensions>>(EMPTY_OBJECT);

  const [isScaled, setIsScaled] = useState(true);

  const {
    data: rawCellTypesByTissueName,
    isLoading: isLoadingCellTypesByTissueName,
  } = useCellTypesByTissueName();

  const [cellTypesByTissueName, setCellTypesByTissueName] =
    useState<CellTypeByTissueName>(EMPTY_OBJECT);

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

  const selectedGeneExpressionSummariesByTissueName = useMemo(() => {
    const result: { [tissueName: string]: GeneExpressionSummary[] } = {};

    for (const tissueName of Object.keys(cellTypesByTissueName)) {
      const tissueGeneExpressionSummaries =
        geneExpressionSummariesByTissueName[tissueName];

      if (!tissueGeneExpressionSummaries) continue;

      result[tissueName] = selectedGenes.map((geneName) => {
        // (thuang): This is needed to ensure the heatmap's gene column
        // is available even if there's no expression data for the column.
        // Otherwise the heatmap columns and column labels won't match up
        // where there's holes in the data.
        const emptyGeneExpressionSummary = {
          cellTypeGeneExpressionSummaries: EMPTY_ARRAY,
          name: geneName,
        };

        return (
          tissueGeneExpressionSummaries[geneName] || emptyGeneExpressionSummary
        );
      });
    }

    return result;
  }, [
    geneExpressionSummariesByTissueName,
    selectedGenes,
    cellTypesByTissueName,
  ]);

  // Listen to delete keyboard press event
  useEffect(() => {
    document.addEventListener("keydown", handleKeyDown);

    return () => {
      document.removeEventListener("keydown", handleKeyDown);
    };

    function handleKeyDown(event: KeyboardEvent): void {
      if (event.code === "Backspace") {
        if (!dispatch) return;

        dispatch(deleteSelectedGenes());
      }
    }
  }, [dispatch]);

  const hasSelectedTissues = (selectedTissues?.length ?? 0) > 0;
  const hasSelectedGenes = selectedGenes.length > 0;

  const shouldShowHeatMap = useMemo(() => {
    return hasSelectedTissues || hasSelectedGenes;
  }, [hasSelectedTissues, hasSelectedGenes]);

  const [isSourceDatasetSidebarOpen, setSourceDatasetSidebarOpen] =
    useState(false);
  const handleSourceDatasetButtonClick = useCallback(() => {
    setSourceDatasetSidebarOpen(!isSourceDatasetSidebarOpen);
  }, [isSourceDatasetSidebarOpen]);

  const [forceOpen, setForceOpen] = useState(false);

  const [downloadStatus, setDownloadStatus] = useState<{
    isLoading: boolean;
  }>({
    isLoading: false,
  });

  const usePrevious = <T,>(value: T): T | undefined => {
    const ref = useRef<T>();
    useEffect(() => {
      ref.current = value;
    });
    return ref.current;
  };
  const prevState = usePrevious({ cellInfoCellType });
  useEffect(() => {
    if (
      prevState?.cellInfoCellType?.cellType.id !== cellInfoCellType?.cellType.id
    ) {
      setForceOpen(!forceOpen); //the value of this boolean isn't actually read downstream, it just checks for uniqueness across renders
    }
  }, [cellInfoCellType, prevState?.cellInfoCellType?.cellType.id, forceOpen]);

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
        disabled={false}
        forceOpen={true}
        wmgSideBar
      >
        <Filters
          isLoading={isLoading}
          availableFilters={availableFilters}
          setAvailableFilters={setAvailableFilters}
          setIsScaled={setIsScaled}
        />
      </SideBar>
      {cellInfoCellType && tissuesByID ? (
        <RightSideBar width={CELL_INFO_SIDEBAR_WIDTH_PX}>
          <CellInfoBar
            generateGeneInfo={generateGeneInfo}
            cellInfoCellType={cellInfoCellType}
            tissueInfo={tissuesByID[cellInfoCellType.tissueID]}
            handleClose={handleCloseRightSideBar}
            title={`${cellInfoCellType.cellType.name}`}
          />

          {
            // Split right sidebar view if fmg AND gene info is populated
            geneInfoGene && (
              <GeneInfoBar
                geneInfoGene={geneInfoGene}
                handleClose={handleCloseGeneInfoSideBar}
                title={geneInfoGene}
              />
            )
          }
        </RightSideBar>
      ) : (
        // Gene info full right sidebar length
        geneInfoGene && (
          <RightSideBar>
            <GeneInfoBar
              geneInfoGene={geneInfoGene}
              handleClose={handleCloseGeneInfoSideBar}
              title={geneInfoGene}
            />
          </RightSideBar>
        )
      )}

      <View id="view" overflow="hidden">
        <Wrapper>
          {isLoading && !shouldShowHeatMap && <Loader />}

          {/* Used for PNG and SVG exports to render message banner to render in output */}
          {downloadStatus.isLoading && (
            <StyledBannerContainer>
              <UnderlyingDataChangeBanner />
            </StyledBannerContainer>
          )}

          <Top id="top-legend">
            <GeneSearchBar className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME} />
            <Legend
              selectedCellTypes={cellTypesByTissueName}
              selectedGenes={selectedGenes}
              selectedTissues={selectedTissues}
              isScaled={isScaled}
              handleRightSidebarButtonClick={handleSourceDatasetButtonClick}
              setDownloadStatus={setDownloadStatus}
              setEchartsRendererMode={setEchartsRendererMode}
              allChartProps={allChartProps}
              availableFilters={availableFilters}
            />
          </Top>

          <GetStarted
            tissueSelected={hasSelectedTissues}
            isLoading={isLoading}
            geneSelected={hasSelectedGenes}
          />

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

          {shouldShowHeatMap ? (
            <HeatMap
              echartsRendererMode={echartsRendererMode}
              cellTypeSortBy={sortBy.cellTypes}
              geneSortBy={sortBy.genes}
              selectedTissues={selectedTissues ?? EMPTY_ARRAY}
              isScaled={isScaled}
              isLoadingAPI={isLoading}
              cellTypes={cellTypesByTissueName}
              genes={selectedGenes}
              selectedGeneExpressionSummariesByTissueName={
                selectedGeneExpressionSummariesByTissueName
              }
              scaledMeanExpressionMax={scaledMeanExpressionMax}
              scaledMeanExpressionMin={scaledMeanExpressionMin}
              selectedOrganismId={selectedOrganismId}
              allChartProps={allChartProps}
              setAllChartProps={setAllChartProps}
            />
          ) : null}
        </Wrapper>
      </View>

      <BottomBanner includeSurveyLink />
    </>
  );
}
