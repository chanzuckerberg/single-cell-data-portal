import cloneDeep from "lodash/cloneDeep";
import { memo, useContext, useMemo, useRef, useState } from "react";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { useResizeObserver } from "src/common/hooks/useResizeObserver";
import {
  generateTermsByKey,
  OntologyTerm,
  usePrimaryFilterDimensions,
} from "src/common/queries/wheresMyGene";
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
import { useSortedCellTypesByTissueName } from "./hooks/useSortedCellTypesByTissueName";
import {
  useSortedGeneNames,
  useTissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap,
} from "./hooks/useSortedGeneNames";
import { useTrackHeatMapLoaded } from "./hooks/useTrackHeatMapLoaded";
import { ChartWrapper, Container, YAxisWrapper } from "./style";
import { getHeatmapWidth, X_AXIS_CHART_HEIGHT_PX, Y_AXIS_CHART_WIDTH_PX } from "./utils";

export interface SelectedGeneExpressionSummariesByTissueName {
    [groupName: string]: {
      [tissueName: string]: GeneExpressionSummary[]
    };
};

interface Props {
  className?: string;
  selectedTissues: string[];
  cellTypes: { [tissue: Tissue]: CellType[] };
  genes: State["selectedGenes"];
  tissuesWithDeletedCellTypes: string[];
  allTissueCellTypes: { [tissue: Tissue]: CellType[] };
  selectedGeneExpressionSummariesByTissueName: SelectedGeneExpressionSummariesByTissueName;
  scaledMeanExpressionMax: number;
  scaledMeanExpressionMin: number;
  isLoadingAPI: boolean;
  isScaled: boolean;
  cellTypeSortBy: SORT_BY;
  geneSortBy: SORT_BY;
  selectedOrganismId: string;
}

export default memo(function HeatMap({
  className,
  selectedTissues,
  cellTypes,
  genes,
  tissuesWithDeletedCellTypes,
  allTissueCellTypes,
  selectedGeneExpressionSummariesByTissueName,
  scaledMeanExpressionMax,
  scaledMeanExpressionMin,
  isLoadingAPI,
  isScaled,
  cellTypeSortBy,
  geneSortBy,
  selectedOrganismId,
}: Props): JSX.Element {
  useTrackHeatMapLoaded({ selectedGenes: genes, selectedTissues });

  const regularGenes = genes.get("") || [];
  // Loading state per tissue
  const [isLoading, setIsLoading] = useState(setInitialIsLoading(cellTypes));
  const chartWrapperRef = useRef<HTMLDivElement>(null);
  const chartWrapperRect = useResizeObserver(chartWrapperRef);

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

  const sortedGeneNamesByGroupName = useSortedGeneNames({
    geneSortBy,
    genes,
    selectedCellTypes: cellTypes,
    tissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap,
  });

  const sortedCellTypesByTissueName = useSortedCellTypesByTissueName({
    cellTypeSortBy,
    genes: regularGenes,
    selectedCellTypes: cellTypes,
    tissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap,
  });

  const geneNameToIndex = useMemo(() => {
    const result: {[groupName: string]: { [key: string]: number }} = {};
    for (const [groupName, sortedGeneNames] of Object.entries(sortedGeneNamesByGroupName)) {
      result[groupName] = {};
      for (const [index, gene] of Object.entries(sortedGeneNames)) {
        result[groupName][gene] = Number(index);
      }
    }
    return result;
  }, [sortedGeneNamesByGroupName]);

  const orderedSelectedGeneExpressionSummariesByTissueName = useMemo(() => {
    const result: SelectedGeneExpressionSummariesByTissueName = {};
    for (const [groupName, tissueNameToGeneExpressionSummary] of Object.entries(
      selectedGeneExpressionSummariesByTissueName
    )) {
      result[groupName] = {};
      for (const [tissueName, geneExpressionSummary] of Object.entries(
        tissueNameToGeneExpressionSummary
      )) {
        // (thuang): sort() mutates the array, so we need to clone it
        result[groupName][tissueName] = cloneDeep(
          geneExpressionSummary.sort((a, b) => {
            if (!a || !b) return -1;

            return geneNameToIndex[groupName][a.name] - geneNameToIndex[groupName][b.name];
          })
        );
      }
    }
    return result;
  }, [selectedGeneExpressionSummariesByTissueName, geneNameToIndex]);

  const geneGroups = Object.entries(sortedGeneNamesByGroupName).map(([_, sortedGeneNames]) => {
    return sortedGeneNames;
  });
  geneGroups.reverse();
  const heatmapOffsets = useMemo(() => {
    const result: number[] = [0];
    for (const [_, sortedGeneNames] of Object.entries(sortedGeneNamesByGroupName).slice().reverse()) {
      result.push(getHeatmapWidth(sortedGeneNames));
      result[result.length - 1] += result[result.length - 2] + 40;
      if (result.length === geneGroups.length) {
        break;
      }
    }
    return result;
  }, [sortedGeneNamesByGroupName, geneGroups]);


  return (
    <Container {...{ className }}>
      {isLoadingAPI || isAnyTissueLoading(isLoading) ? <Loader /> : null}
      <div
        style={{
          display: "flex",
          backgroundColor: "white",
          flexDirection: "row",
          position: "sticky",
          top: 0,
          width: "100%",
          zIndex: 2
        }}
      >
        <div style={{width: Y_AXIS_CHART_WIDTH_PX, height: X_AXIS_CHART_HEIGHT_PX}}/>
        <CellCountLabel>Cell Count</CellCountLabel>
          {geneGroups.map(
            (sortedGeneNames, index) => {
              return (
                <XAxisChart 
                  key={`${index}-x-axis-chart`}
                  geneNames={sortedGeneNames}
                  leftOffset={Y_AXIS_CHART_WIDTH_PX+heatmapOffsets[index]}
                />
              );
            }
          )}
      </div>  

      <YAxisWrapper
        height={(chartWrapperRect?.height || 0) - X_AXIS_CHART_HEIGHT_PX}
      >
        {selectedTissues.map((tissue) => {
          const tissueCellTypes = getTissueCellTypes({
            cellTypeSortBy,
            cellTypes,
            sortedCellTypesByTissueName,
            tissue,
          });
          
          return (
            <YAxisChart
              key={tissue}
              tissue={tissue}
              tissueID={tissuesByName[tissue].id}
              cellTypes={tissueCellTypes}
              hasDeletedCellTypes={tissuesWithDeletedCellTypes.includes(tissue)}
              availableCellTypes={allTissueCellTypes[tissue]}
              generateMarkerGenes={generateMarkerGenes}
              selectedOrganismId={selectedOrganismId}
            />
          );
        })}
      </YAxisWrapper>
      <ChartWrapper ref={chartWrapperRef}>
        {selectedTissues.map((tissue) => {
          const tissueCellTypes = getTissueCellTypes({
            cellTypeSortBy,
            cellTypes,
            sortedCellTypesByTissueName,
            tissue,
          });
          const els: JSX.Element[] = [];
          Object.entries(orderedSelectedGeneExpressionSummariesByTissueName).forEach((
            [groupName, orderedSelectedGeneExpressionSummaries], index
          )=>{
            els.push(
              <Chart
                isScaled={isScaled}
                key={`${tissue}-${groupName}`}
                tissue={tissue}
                cellTypes={tissueCellTypes}
                selectedGeneData={
                  orderedSelectedGeneExpressionSummaries[tissue]
                }
                setIsLoading={setIsLoading}
                scaledMeanExpressionMax={scaledMeanExpressionMax}
                scaledMeanExpressionMin={scaledMeanExpressionMin}
                leftOffset={heatmapOffsets[index]}
              />
            );
          });
          els.reverse();
          return (
            <div
              key={`${tissue}-chart`}
              style={{ 
                display: "flex",
                flexDirection: "row",
                position: "relative",
                left: Y_AXIS_CHART_WIDTH_PX,
                columnGap: "40px"
              }}
            >   
              {els}               
            </div>   
          );
        })}
      </ChartWrapper>
    </Container>
  );
});

function getTissueCellTypes({
  cellTypes,
  sortedCellTypesByTissueName,
  tissue,
  cellTypeSortBy,
}: {
  cellTypes: { [tissue: Tissue]: CellType[] };
  sortedCellTypesByTissueName: { [tissue: string]: CellType[] };
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
