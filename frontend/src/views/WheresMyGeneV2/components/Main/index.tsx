import { DrawerSize } from "@blueprintjs/core";
import Head from "next/head";
import {
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
  OntologyTerm,
  useCellTypesByTissueName,
  useGeneExpressionSummariesByTissueName,
  usePrimaryFilterDimensions,
} from "src/common/queries/wheresMyGene";
import SideBar from "src/components/common/SideBar";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import {
  addGeneInfoGene,
  clearGeneInfoGene,
  closeRightSidebar,
  deleteSelectedGenes,
} from "src/views/WheresMyGene/common/store/actions";
import {
  GeneExpressionSummary,
  Tissue,
  ChartProps,
} from "src/views/WheresMyGene/common/types";
import CellInfoSideBar from "src/views/WheresMyGene/components/CellInfoSideBar";
import Filters from "src/views/WheresMyGene/components/Filters";
import GeneInfoSideBar from "src/views/WheresMyGene/components/GeneInfoSideBar";
import GeneSearchBar from "src/views/WheresMyGeneV2/components/GeneSearchBar";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "src/views/WheresMyGene/components/GeneSearchBar/components/SaveExport";
import { UnderlyingDataChangeBanner } from "src/views/WheresMyGene/components/GeneSearchBar/components/SaveExport/ExportBanner";

import InfoPanel from "src/views/WheresMyGene/components/InfoPanel";
import Legend from "src/views/WheresMyGene/components/InfoPanel/components/Legend";
import Loader from "src/views/WheresMyGene/components/Loader";
import {
  StyledBannerContainer,
  StyledSidebarDrawer,
} from "src/views/WheresMyGene/components/Main/style";
import RightSideBar from "src/views/WheresMyGene/components/RightSideBar";
import ScreenTint from "src/views/WheresMyGene/components/ScreenTint";
import {
  SideBarPositioner,
  SideBarWrapper,
  Top,
  Wrapper,
} from "src/views/WheresMyGene/style";
import { View } from "src/views/globalStyle";
import HeatMap from "../HeatMap";
import BottomBanner from "src/components/BottomBanner";

export const INFO_PANEL_WIDTH_PX = 320;

export default function WheresMyGene(): JSX.Element {
  const state = useContext(StateContext);
  const dispatch = useContext(DispatchContext);

  const { selectedGenes, sortBy, geneInfoGene, cellInfoCellType } = state;

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

  // This is set in HeatMap and the value is used to determine spacing in SVG export
  const [expandedTissues, setExpandedTissues] = useState<Set<Tissue>>(
    new Set()
  );

  //(seve): These useEffects are deceptively simple.
  // Their purpose is to avoid updating the state with null/empty values while we're waiting for the api to return data.

  const {
    data: rawCellTypesByTissueName,
    isLoading: isLoadingCellTypesByTissueName,
  } = useCellTypesByTissueName(2);

  const [cellTypesByTissueName, setCellTypesByTissueName] =
    useState<CellTypeByTissueName>(EMPTY_OBJECT);

  // This is needed
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

  const sortedTissues = useMemo(() => {
    return Object.keys(tissuesByName);
  }, [tissuesByName]);

  return (
    <>
      <Head>
        <title>CELL&times;GENE | Gene Expression</title>
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
        <RightSideBar>
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
        </RightSideBar>
      ) : (
        // Gene info full right sidebar length
        geneInfoGene && (
          <RightSideBar>
            <GeneInfoSideBar
              geneInfoGene={geneInfoGene}
              handleClose={handleCloseGeneInfoSideBar}
              title={`${geneInfoGene}`}
            />
          </RightSideBar>
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
            <GeneSearchBar className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME} />
            <Legend
              selectedCellTypes={cellTypesByTissueName}
              selectedGenes={selectedGenes}
              isScaled={isScaled}
              handleRightSidebarButtonClick={handleSourceDatasetButtonClick}
              setDownloadStatus={setDownloadStatus}
              setEchartsRendererMode={setEchartsRendererMode}
              allChartProps={allChartProps}
              availableFilters={availableFilters}
              tissues={sortedTissues}
              expandedTissues={expandedTissues}
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
            selectedOrganismId={selectedOrganismId}
            allChartProps={allChartProps}
            setAllChartProps={setAllChartProps}
            setTissuesByName={setTissuesByName}
            expandedTissues={expandedTissues}
            setExpandedTissues={setExpandedTissues}
          />
        </Wrapper>
      </View>

      <BottomBanner includeSurveyLink />
    </>
  );
}
