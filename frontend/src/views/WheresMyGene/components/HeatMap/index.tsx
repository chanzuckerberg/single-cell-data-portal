import cloneDeep from "lodash/cloneDeep";
import {
  Dispatch,
  memo,
  SetStateAction,
  useContext,
  useMemo,
  useRef,
  useState,
} from "react";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import {
  CellTypeRow,
  generateTermsByKey,
  OntologyTerm,
  usePrimaryFilterDimensions,
} from "src/common/queries/wheresMyGene";
import {
  HEATMAP_CONTAINER_ID,
  X_AXIS_CHART_HEIGHT_PX,
} from "../../../WheresMyGeneV2/common/constants";
import {
  DispatchContext,
  State,
  StateContext,
} from "../../../WheresMyGeneV2/common/store";
import { addCellInfoCellType } from "../../../WheresMyGeneV2/common/store/actions";
import {
  CellType,
  ChartProps,
  GeneExpressionSummary,
  SORT_BY,
  Tissue,
} from "../../../WheresMyGeneV2/common/types";
import Loader from "../Loader";
import Chart from "../../../WheresMyGeneV2/components/HeatMap/components/Chart";
import XAxisChart from "../../../WheresMyGeneV2/components/HeatMap/components/XAxisChart";
import { CellCountLabel } from "../../../WheresMyGeneV2/components/HeatMap/components/XAxisChart/style";
import YAxisChart from "./components/YAxisChart";
import { useSortedCellTypesByTissueName } from "./hooks/useSortedCellTypesByTissueName";
import {
  useSortedGeneNames,
  useTissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap,
} from "./hooks/useSortedGeneNames";
import { useTrackHeatMapLoaded } from "./hooks/useTrackHeatMapLoaded";
import {
  ChartWrapper,
  Container,
  ContainerWrapper,
  TopLeftCornerMask,
  XAxisMask,
  XAxisWrapper,
  YAxisWrapper,
} from "./style";
import { HEAT_MAP_BASE_CELL_PX, hyphenize } from "./utils";

interface Props {
  className?: string;
  selectedTissues: string[];
  cellTypes: { [tissue: Tissue]: CellTypeRow[] };
  genes: State["selectedGenes"];
  selectedGeneExpressionSummariesByTissueName: {
    [tissueName: string]: GeneExpressionSummary[];
  };
  scaledMeanExpressionMax: number;
  scaledMeanExpressionMin: number;
  isLoadingAPI: boolean;
  isScaled: boolean;
  cellTypeSortBy: SORT_BY;
  geneSortBy: SORT_BY;
  selectedOrganismId: string;
  echartsRendererMode: "svg" | "canvas";
  setAllChartProps: Dispatch<
    SetStateAction<{
      [tissue: string]: ChartProps;
    }>
  >;
  allChartProps: { [tissue: string]: ChartProps };
}

export default memo(function HeatMap({
  className,
  selectedTissues,
  cellTypes,
  genes,
  selectedGeneExpressionSummariesByTissueName,
  scaledMeanExpressionMax,
  scaledMeanExpressionMin,
  isLoadingAPI,
  isScaled,
  cellTypeSortBy,
  geneSortBy,
  selectedOrganismId,
  echartsRendererMode,
  allChartProps,
  setAllChartProps,
}: Props): JSX.Element {
  useTrackHeatMapLoaded({ selectedGenes: genes, selectedTissues });

  const { xAxisHeight } = useContext(StateContext);

  // Loading state per tissue
  const [isLoading, setIsLoading] = useState(setInitialIsLoading(cellTypes));
  const chartWrapperRef = useRef<HTMLDivElement>(null);
  const dispatch = useContext(DispatchContext);

  const { data } = usePrimaryFilterDimensions();

  // Get tissueName to ID map for use in find marker genes
  const tissuesByName = useMemo(() => {
    let result: { [name: string]: OntologyTerm } = {};

    if (!data) return result;

    const { tissues } = data;

    result = generateTermsByKey(tissues, "name");

    return result;
  }, [data]);

  const generateMarkerGenes = (cellType: CellType, tissueID: string) => {
    if (!dispatch) return;
    dispatch(addCellInfoCellType({ cellType, tissueID }));
  };

  const tissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap =
    useTissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap(
      selectedGeneExpressionSummariesByTissueName
    );

  const sortedGeneNames = useSortedGeneNames({
    geneSortBy,
    genes,
    selectedCellTypes: cellTypes,
    tissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap,
  });

  const sortedCellTypesByTissueName = useSortedCellTypesByTissueName({
    cellTypeSortBy,
    genes,
    selectedCellTypes: cellTypes,
    tissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap,
  });

  const geneNameToIndex = useMemo(() => {
    const result: { [key: string]: number } = {};

    for (const [index, gene] of Object.entries(sortedGeneNames)) {
      result[gene] = Number(index);
    }

    return result;
  }, [sortedGeneNames]);

  const orderedSelectedGeneExpressionSummariesByTissueName = useMemo(() => {
    const result: { [tissueName: string]: GeneExpressionSummary[] } = {};

    for (const [tissueName, geneExpressionSummary] of Object.entries(
      selectedGeneExpressionSummariesByTissueName
    )) {
      // (thuang): sort() mutates the array, so we need to clone it
      result[tissueName] = cloneDeep(
        geneExpressionSummary.sort((a, b) => {
          if (!a || !b) return -1;

          return geneNameToIndex[a.name] - geneNameToIndex[b.name];
        })
      );
    }

    return result;
  }, [selectedGeneExpressionSummariesByTissueName, geneNameToIndex]);

  return (
    <ContainerWrapper>
      <TopLeftCornerMask height={xAxisHeight}>
        <CellCountLabel>Cell Count</CellCountLabel>
      </TopLeftCornerMask>
      <Container {...{ className }} id={HEATMAP_CONTAINER_ID}>
        {isLoadingAPI || isAnyTissueLoading(isLoading) ? <Loader /> : null}
        <XAxisWrapper id="x-axis-wrapper">
          <XAxisMask data-testid="x-axis-mask" height={xAxisHeight} />
          <XAxisChart
            geneNames={sortedGeneNames}
            // (thuang): Dummy value here, since we're deleting this file
            sidebarWidth={0}
          />
        </XAxisWrapper>
        <YAxisWrapper top={xAxisHeight}>
          {selectedTissues.map((tissue) => {
            const tissueCellTypes = getTissueCellTypes({
              cellTypeSortBy,
              cellTypes,
              sortedCellTypesByTissueName,
              tissue,
            });
            return tissueCellTypes.length ? (
              <div id={`y-axis-${hyphenize(tissue)}`}>
                <YAxisChart
                  key={tissue}
                  tissue={tissue}
                  tissueID={tissuesByName[tissue].id}
                  cellTypes={tissueCellTypes}
                  generateMarkerGenes={generateMarkerGenes}
                  selectedOrganismId={selectedOrganismId}
                />
              </div>
            ) : null;
          })}
        </YAxisWrapper>
        <ChartWrapper ref={chartWrapperRef} top={xAxisHeight}>
          {selectedTissues.map((tissue) => {
            const tissueCellTypes = getTissueCellTypes({
              cellTypeSortBy,
              cellTypes,
              sortedCellTypesByTissueName,
              tissue,
            });

            const selectedGeneData =
              orderedSelectedGeneExpressionSummariesByTissueName[tissue];

            /**
             * (thuang): If there is no selected gene data, we don't want to render
             * the chart, because it will cause the chart to render with 0 width,
             * which is an error for echarts
             */

            if (!selectedGeneData?.length) {
              const height = tissueCellTypes.length * HEAT_MAP_BASE_CELL_PX;
              return (
                <div
                  key={`y-axis-${hyphenize(tissue)}`}
                  style={{
                    // If the height is 0, then set the heatmap height as 0. Else render the height normally
                    height: `${height && height + X_AXIS_CHART_HEIGHT_PX}px`,
                  }}
                />
              );
            }

            return (
              <Chart
                isScaled={isScaled}
                /**
                 * (thuang): We use `key` to force re-render the HeatMap component
                 * when the renderer mode changes, so echarts can create new instances
                 */
                key={`${tissue}-${echartsRendererMode}`}
                tissue={tissue}
                cellTypes={tissueCellTypes}
                selectedGeneData={
                  orderedSelectedGeneExpressionSummariesByTissueName[tissue]
                }
                setIsLoading={setIsLoading}
                scaledMeanExpressionMax={scaledMeanExpressionMax}
                scaledMeanExpressionMin={scaledMeanExpressionMin}
                echartsRendererMode={echartsRendererMode}
                setAllChartProps={setAllChartProps}
                chartProps={allChartProps[tissue]}
                maxExpression={6.0}
              />
            );
          })}
        </ChartWrapper>
      </Container>
    </ContainerWrapper>
  );
});

function getTissueCellTypes({
  cellTypes,
  sortedCellTypesByTissueName,
  tissue,
  cellTypeSortBy,
}: {
  cellTypes: { [tissue: Tissue]: CellTypeRow[] };
  sortedCellTypesByTissueName: { [tissue: string]: CellTypeRow[] };
  tissue: Tissue;
  cellTypeSortBy: SORT_BY;
}) {
  const tissueCellTypes = cellTypes[tissue];
  const sortedTissueCellTypes = sortedCellTypesByTissueName[tissue];
  return (
    (cellTypeSortBy === SORT_BY.CELL_ONTOLOGY
      ? tissueCellTypes
      : sortedTissueCellTypes) || EMPTY_ARRAY
  );
}

function isAnyTissueLoading(isLoading: { [tissue: Tissue]: boolean }) {
  return Object.values(isLoading).some((isLoading) => isLoading);
}

function setInitialIsLoading(cellTypes: Props["cellTypes"]) {
  return Object.keys(cellTypes).reduce((isLoading, tissue) => {
    return { ...isLoading, [tissue]: false };
  }, {});
}
