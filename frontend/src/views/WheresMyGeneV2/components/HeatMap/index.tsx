import { Button } from "czifui";
import cloneDeep from "lodash/cloneDeep";
import {
  Dispatch,
  memo,
  SetStateAction,
  useContext,
  useEffect,
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
} from "src/common/queries/wheresMyGeneV2";
import { HEATMAP_CONTAINER_ID } from "../../common/constants";
import { DispatchContext, State } from "../../common/store";
import { addCellInfoCellType } from "../../common/store/actions";
import {
  CellType,
  GeneExpressionSummary,
  SORT_BY,
  Tissue,
} from "../../common/types";
import Loader from "../Loader";
import Chart from "./components/Chart";
import XAxisChart from "./components/XAxisChart";
import { CellCountLabel } from "./components/XAxisChart/style";
import YAxisChart from "./components/YAxisChart";
import { ChartProps } from "./hooks/common/types";
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

interface Props {
  className?: string;
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
  useTrackHeatMapLoaded({ selectedGenes: genes });

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

  // would this be more performant as a bitmap of booleans?
  const initialBitMap = useMemo(() => {
    return Object.entries(sortedCellTypesByTissueName).reduce(
      (acc, [tissue, cellTypes]) => {
        acc[tissue] = {};
        cellTypes.forEach((cellType) => {
          acc[tissue][cellType.id] = false;
        });
        acc[tissue]["HEADER_ROW"] = true;
        return acc;
      },
      {} as { [tissue: string]: { [cellType: string]: boolean } }
    );
  }, [sortedCellTypesByTissueName]);
  // (seve): This state should be moved to the store
  // This is a 2 dimensional array of booleans that represents which rows belonging to each tissue displayed
  const [displayedRowBitMap, setDisplayedRowBitMap] = useState(initialBitMap);

  useEffect(() => {
    setDisplayedRowBitMap(initialBitMap);
  }, [initialBitMap]);

  const handleExpand = () => {
    //expand the first tissue
    displayedRowBitMap["lung"] = Object.fromEntries(
      Object.entries(displayedRowBitMap["lung"]).map(([key]) => {
        return [key, true];
      })
    );
    setDisplayedRowBitMap(displayedRowBitMap);
  };

  return (
    <>
      <Button onClick={handleExpand}>EXPAND</Button>

      <ContainerWrapper>
        <TopLeftCornerMask>
          <CellCountLabel>Cell Count</CellCountLabel>
        </TopLeftCornerMask>
        <Container {...{ className }} id={HEATMAP_CONTAINER_ID}>
          {isLoadingAPI || isAnyTissueLoading(isLoading) ? <Loader /> : null}
          <XAxisWrapper id="x-axis-wrapper">
            <XAxisMask data-testid="x-axis-mask" />
            <XAxisChart geneNames={sortedGeneNames} />
          </XAxisWrapper>
          <YAxisWrapper>
            {Object.values(tissuesByName).map((tissue: OntologyTerm) => {
              const tissueCellTypes = getTissueCellTypes({
                cellTypeSortBy,
                cellTypes,
                sortedCellTypesByTissueName,
                tissue: tissue.name,
                displayedRowBitMap,
              });
              return (
                <YAxisChart
                  key={tissue.name}
                  tissue={tissue.name}
                  tissueID={tissue.id}
                  cellTypes={tissueCellTypes}
                  generateMarkerGenes={generateMarkerGenes}
                  selectedOrganismId={selectedOrganismId}
                />
              );
            })}
          </YAxisWrapper>
          <ChartWrapper ref={chartWrapperRef}>
            {Object.values(tissuesByName).map((tissue: OntologyTerm) => {
              const tissueCellTypes = getTissueCellTypes({
                cellTypeSortBy,
                cellTypes,
                sortedCellTypesByTissueName,
                tissue: tissue.name,
                displayedRowBitMap,
              });

              const selectedGeneData =
                orderedSelectedGeneExpressionSummariesByTissueName[tissue.name];

              /**
               * (thuang): If there is no selected gene data, we don't want to render
               * the chart, because it will cause the chart to render with 0 width,
               * which is an error for echarts
               */
              if (!selectedGeneData?.length) return null;

              return (
                <Chart
                  isScaled={isScaled}
                  /**
                   * (thuang): We use `key` to force re-render the HeatMap component
                   * when the renderer mode changes, so echarts can create new instances
                   */
                  key={`${tissue.name}-${echartsRendererMode}`}
                  tissue={tissue.name}
                  cellTypes={tissueCellTypes}
                  selectedGeneData={
                    orderedSelectedGeneExpressionSummariesByTissueName[
                      tissue.name
                    ]
                  }
                  setIsLoading={setIsLoading}
                  scaledMeanExpressionMax={scaledMeanExpressionMax}
                  scaledMeanExpressionMin={scaledMeanExpressionMin}
                  echartsRendererMode={echartsRendererMode}
                  setAllChartProps={setAllChartProps}
                  chartProps={allChartProps[tissue.name]}
                />
              );
            })}
          </ChartWrapper>
        </Container>
      </ContainerWrapper>
    </>
  );
});

function getTissueCellTypes({
  cellTypes,
  sortedCellTypesByTissueName,
  tissue,
  cellTypeSortBy,
  displayedRowBitMap,
}: {
  cellTypes: { [tissue: Tissue]: CellTypeRow[] };
  sortedCellTypesByTissueName: { [tissue: string]: CellTypeRow[] };
  tissue: Tissue;
  cellTypeSortBy: SORT_BY;
  displayedRowBitMap: { [tissue: string]: { [cellType: string]: boolean } };
}) {
  const tissueCellTypes = cellTypes[tissue];
  const sortedTissueCellTypes = sortedCellTypesByTissueName[tissue];
  const ret =
    (cellTypeSortBy === SORT_BY.CELL_ONTOLOGY
      ? tissueCellTypes
      : sortedTissueCellTypes) || EMPTY_ARRAY;
  if (ret[0]?.id !== "HEADER_ROW")
    ret.unshift({ ...ret[0], id: "HEADER_ROW", name: tissue });

  return ret.filter((cellType) => displayedRowBitMap[tissue]?.[cellType.id]);
}

function isAnyTissueLoading(isLoading: { [tissue: Tissue]: boolean }) {
  return Object.values(isLoading).some((isLoading) => isLoading);
}

function setInitialIsLoading(cellTypes: Props["cellTypes"]) {
  return Object.keys(cellTypes).reduce((isLoading, tissue) => {
    return { ...isLoading, [tissue]: false };
  }, {});
}
