import { Button } from "@czi-sds/components";
import cloneDeep from "lodash/cloneDeep";
import {
  Dispatch,
  memo,
  SetStateAction,
  useCallback,
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
import { memoize } from "lodash";

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

  const initialDisplayedCellTypes = useMemo(() => {
    return Object.entries(sortedCellTypesByTissueName).reduce(
      (acc, [tissue]) => {
        acc.add(tissue + tissue);
        return acc;
      },
      new Set<string>()
    );
  }, [sortedCellTypesByTissueName]);
  // (seve): This state should be moved to the store
  const [displayedCellTypes, setDisplayedCellTypes] = useState<Set<string>>(
    initialDisplayedCellTypes
  );

  // set of tissue names that are visible and set of cell types that are visible
  // presence is visible

  useEffect(() => {
    setDisplayedCellTypes(initialDisplayedCellTypes);
  }, [initialDisplayedCellTypes]);

  const [expandedTissues, setExpandedTissues] = useState<Array<Tissue>>([]);

  const handleExpand = useCallback(() => {
    const newDisplayedCellTypes = new Set<string>(displayedCellTypes);

    if (expandedTissues.includes("lung")) {
      setExpandedTissues([]);
      setDisplayedCellTypes(initialDisplayedCellTypes);
      console.log("DEFAULT");

      return;
    }
    //expand the first tissue

    sortedCellTypesByTissueName["lung"].forEach((cellType) => {
      newDisplayedCellTypes.add("lung" + cellType.name);
    });
    setExpandedTissues(["lung"]);
    setDisplayedCellTypes(newDisplayedCellTypes);
  }, [
    displayedCellTypes,
    sortedCellTypesByTissueName,
    expandedTissues,
    initialDisplayedCellTypes,
  ]);

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
              const tissueCellTypes = memoizedGetTissueCellTypes({
                cellTypeSortBy,
                cellTypes,
                sortedCellTypesByTissueName,
                tissue: tissue.name,
                displayedCellTypes: displayedCellTypes,
              });

              return (
                tissueCellTypes.length > 0 && (
                  <YAxisChart
                    key={tissue.name}
                    tissue={tissue.name}
                    tissueID={tissue.id}
                    cellTypes={tissueCellTypes}
                    generateMarkerGenes={generateMarkerGenes}
                    selectedOrganismId={selectedOrganismId}
                  />
                )
              );
            })}
          </YAxisWrapper>
          <ChartWrapper ref={chartWrapperRef}>
            {Object.values(tissuesByName).map((tissue: OntologyTerm) => {
              const tissueCellTypes = memoizedGetTissueCellTypes({
                cellTypeSortBy,
                cellTypes,
                sortedCellTypesByTissueName,
                tissue: tissue.name,
                displayedCellTypes: displayedCellTypes,
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
  displayedCellTypes,
}: {
  cellTypes: { [tissue: Tissue]: CellTypeRow[] };
  sortedCellTypesByTissueName: { [tissue: string]: CellTypeRow[] };
  tissue: Tissue;
  cellTypeSortBy: SORT_BY;
  displayedCellTypes: Set<string>;
}) {
  const tissueCellTypes = cellTypes[tissue];
  if (!tissueCellTypes || tissueCellTypes.length === 0) return [];
  const sortedTissueCellTypes = sortedCellTypesByTissueName[tissue];
  let ret =
    (cellTypeSortBy === SORT_BY.CELL_ONTOLOGY
      ? tissueCellTypes
      : sortedTissueCellTypes) || EMPTY_ARRAY;

  if (ret[ret.length - 1].name !== tissue)
    ret.push({ ...ret[0], name: tissue });

  ret = ret.filter((cellType) =>
    displayedCellTypes.has(tissue + cellType.name)
  );
  return ret;
}

const memoizedGetTissueCellTypes = memoize(
  getTissueCellTypes,
  ({
    cellTypes,
    sortedCellTypesByTissueName,
    tissue,
    cellTypeSortBy,
    displayedCellTypes,
  }) => {
    return `${tissue}-${cellTypeSortBy}-${[...displayedCellTypes]?.join(
      ""
    )}-${cellTypes[tissue]?.join("")}-${sortedCellTypesByTissueName[
      tissue
    ]?.join("")}`;
  }
);

function isAnyTissueLoading(isLoading: { [tissue: Tissue]: boolean }) {
  return Object.values(isLoading).some((isLoading) => isLoading);
}

function setInitialIsLoading(cellTypes: Props["cellTypes"]) {
  return Object.keys(cellTypes).reduce((isLoading, tissue) => {
    return { ...isLoading, [tissue]: false };
  }, {});
}
