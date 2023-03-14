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
  OntologyTerm,
  useCellTypesByTissueName,
  useGeneExpressionSummariesByTissueName,
  usePrimaryFilterDimensions,
} from "src/common/queries/wheresMyGene";
import SideBar from "src/components/common/SideBar";
import {
  GeneSideBarOpenButtonWrapper,
  Position,
} from "src/components/common/SideBar/style";
import { View } from "../../../globalStyle";
import { DispatchContext, StateContext } from "../../common/store";
import { deleteSelectedGenes } from "../../common/store/actions";
import { GeneExpressionSummary } from "../../common/types";
import { SideBarPositioner, SideBarWrapper, Top, Wrapper } from "../../style";
import Beta from "../Beta";
import CellInfoBar from "../CellInfoSideBar";
import { CELL_INFO_SIDEBAR_WIDTH_PX } from "../CellInfoSideBar/style";
import Filters from "../Filters";
import GeneSearchBar from "../GeneSearchBar";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../GeneSearchBar/components/SaveImage";
import GetStarted from "../GetStarted";
import HeatMap from "../HeatMap";
import { ChartProps } from "../HeatMap/hooks/common/types";
import InfoPanel from "../InfoPanel";
import ColorScale from "../InfoPanel/components/ColorScale";
import Legend from "../InfoPanel/components/Legend";
import Loader from "../Loader";
import ScreenTint from "../ScreenTint";
import { BetaWrapper, SideBarLabel, StyledSidebarDrawer } from "./style";

export const INFO_PANEL_WIDTH_PX = 320;

export default function WheresMyGene(): JSX.Element {
  const state = useContext(StateContext);
  const dispatch = useContext(DispatchContext);

  const { selectedGenes, selectedTissues, sortBy, cellInfoCellType } = state;

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

  const [availableOrganisms, setAvailableOrganisms] =
    useState<OntologyTerm[]>(EMPTY_ARRAY);

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
        (cellType) => cellType.id
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
                cellTypeGeneExpressionSummary.id
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

  const hasSelectedTissues = selectedTissues.length > 0;
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
    blur?: boolean;
  }>({
    isLoading: false,
    blur: false,
  });

  const usePrevious = <T,>(value: T): T | undefined => {
    const ref = useRef<T>();
    useEffect(() => {
      ref.current = value;
    });
    return ref.current;
  };
  const prevAmount = usePrevious({ cellInfoCellType });
  useEffect(() => {
    if (
      prevAmount?.cellInfoCellType?.cellType.id !==
      cellInfoCellType?.cellType.id
    ) {
      setForceOpen(!forceOpen); //the value of this boolean isn't actually read downstream, it just checks for uniqueness across renders
    }
  }, [cellInfoCellType, prevAmount?.cellInfoCellType?.cellType.id, forceOpen]);

  const [echartsRendererMode, setEchartsRendererMode] = useState<
    "canvas" | "svg"
  >("canvas");

  return (
    <>
      <Head>
        <title>CELL&times;GENE | Gene Expression</title>
      </Head>

      <SideBar
        label={<SideBarLabel>Filters</SideBarLabel>}
        SideBarWrapperComponent={SideBarWrapper}
        SideBarPositionerComponent={SideBarPositioner}
        testId="filters-panel"
        disabled={false}
        forceToggle={true}
        wmgSideBar
      >
        <Filters
          isLoading={isLoading}
          availableFilters={availableFilters}
          setAvailableFilters={setAvailableFilters}
          setAvailableOrganisms={setAvailableOrganisms}
        />

        <ColorScale setIsScaled={setIsScaled} />
      </SideBar>
      {cellInfoCellType && tissuesByID && (
        <SideBar
          label={`${cellInfoCellType.cellType.name}`}
          SideBarWrapperComponent={SideBarWrapper}
          SideBarPositionerComponent={SideBarPositioner}
          SideBarOpenButtonWrapperComponent={GeneSideBarOpenButtonWrapper}
          position={Position.RIGHT}
          testId="cell-type-details-panel"
          disabled={false}
          forceToggle={forceOpen}
          wmgSideBar
          width={CELL_INFO_SIDEBAR_WIDTH_PX}
          truncatedLabel={`${tissuesByID[cellInfoCellType.tissueID].name} - ${
            cellInfoCellType.cellType.name
          }`}
        >
          <CellInfoBar
            cellInfoCellType={cellInfoCellType}
            tissueInfo={tissuesByID[cellInfoCellType.tissueID]}
          />
        </SideBar>
      )}

      <View id="view" overflow="hidden">
        <Wrapper>
          {isLoading && !shouldShowHeatMap && <Loader />}

          <Top>
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
              availableOrganisms={availableOrganisms}
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
              selectedTissues={selectedTissues}
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
        <BetaWrapper>
          <Beta className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME} />
        </BetaWrapper>
      </View>
    </>
  );
}
