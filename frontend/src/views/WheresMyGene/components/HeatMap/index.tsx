import { memo, useRef, useState } from "react";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { useResizeObserver } from "src/common/hooks/useResizeObserver";
import { State } from "../../common/store";
import { CellType, GeneExpressionSummary, Tissue } from "../../common/types";
import Loader from "../Loader";
import Chart from "./components/Chart";
import XAxisChart from "./components/XAxisChart";
import YAxisChart from "./components/YAxisChart";
import { useTrackHeatMapLoaded } from "./hooks/useTrackHeatMapLoaded";
import { ChartWrapper, Container, YAxisWrapper } from "./style";
import { X_AXIS_CHART_HEIGHT_PX } from "./utils";

interface Props {
  selectedTissues: string[];
  cellTypes: { [tissue: Tissue]: CellType[] };
  genes: State["selectedGenes"];
  tissuesWithDeletedCellTypes: string[];
  allTissueCellTypes: { [tissue: Tissue]: CellType[] };
  selectedGeneExpressionSummariesByTissueName: {
    [tissueName: string]: GeneExpressionSummary[];
  };
  scaledMeanExpressionMax: number;
  scaledMeanExpressionMin: number;
  isLoadingAPI: boolean;
  isScaled: boolean;
}

export default memo(function HeatMap({
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
}: Props): JSX.Element {
  useTrackHeatMapLoaded({ selectedTissues, selectedGenes: genes });

  // Loading state per tissue
  const [isLoading, setIsLoading] = useState(setInitialIsLoading(cellTypes));
  const chartWrapperRef = useRef<HTMLDivElement>(null);
  const chartWrapperRect = useResizeObserver(chartWrapperRef);

  return (
    <Container>
      {isLoadingAPI || isAnyTissueLoading(isLoading) ? <Loader /> : null}

      <XAxisChart geneNames={genes} />

      <YAxisWrapper
        height={(chartWrapperRect?.height || 0) - X_AXIS_CHART_HEIGHT_PX}
      >
        {selectedTissues.map((tissue) => {
          const tissueCellTypes = cellTypes[tissue];

          return (
            <YAxisChart
              key={tissue}
              tissue={tissue}
              cellTypes={tissueCellTypes}
              hasDeletedCellTypes={tissuesWithDeletedCellTypes.includes(tissue)}
              availableCellTypes={allTissueCellTypes[tissue]}
            />
          );
        })}
      </YAxisWrapper>
      <ChartWrapper ref={chartWrapperRef}>
        {selectedTissues.map((tissue) => {
          const tissueCellTypes = cellTypes[tissue] || EMPTY_ARRAY;

          return (
            <Chart
              isScaled={isScaled}
              key={tissue}
              tissue={tissue}
              cellTypes={tissueCellTypes}
              selectedGeneData={
                selectedGeneExpressionSummariesByTissueName[tissue]
              }
              setIsLoading={setIsLoading}
              scaledMeanExpressionMax={scaledMeanExpressionMax}
              scaledMeanExpressionMin={scaledMeanExpressionMin}
            />
          );
        })}
      </ChartWrapper>
    </Container>
  );
});

function isAnyTissueLoading(isLoading: { [tissue: Tissue]: boolean }) {
  return Object.values(isLoading).some((isLoading) => isLoading);
}

function setInitialIsLoading(cellTypes: Props["cellTypes"]) {
  return Object.keys(cellTypes).reduce((isLoading, tissue) => {
    return { ...isLoading, [tissue]: false };
  }, {});
}
