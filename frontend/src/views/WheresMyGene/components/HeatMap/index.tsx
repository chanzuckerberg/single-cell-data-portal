import { Intent, Spinner } from "@blueprintjs/core";
import { memo, useMemo, useState } from "react";
import { State } from "../../common/store";
import { CellType, GeneExpressionSummary, Tissue } from "../../common/types";
import Chart from "./components/Chart";
import XAxisChart from "./components/XAxisChart";
import YAxisChart from "./components/YAxisChart";
import { ChartWrapper, Container, Loader, YAxisWrapper } from "./style";
import { getHeatmapHeight, X_AXIS_CHART_HEIGHT_PX } from "./utils";

interface Props {
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
  // Loading state per tissue
  const [isLoading, setIsLoading] = useState(setInitialIsLoading(cellTypes));

  const yAxisWrapperHeight = useMemo(() => {
    const yAxisChartHeight = Object.values(cellTypes).reduce(
      (height, cellTypeSummaries) => {
        return height + getHeatmapHeight(cellTypeSummaries);
      },
      0
    );

    const tissueCount = Object.keys(cellTypes).length;

    return yAxisChartHeight * tissueCount + X_AXIS_CHART_HEIGHT_PX;
  }, [cellTypes]);

  return (
    <Container>
      {isLoadingAPI || isAnyTissueLoading(isLoading) ? (
        <Loader>
          <Spinner intent={Intent.PRIMARY} size={20} />
          Loading...
        </Loader>
      ) : null}

      <XAxisChart geneNames={genes} />

      <YAxisWrapper height={yAxisWrapperHeight}>
        {Object.entries(cellTypes).map(([tissue, tissueCellTypes]) => {
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
      <ChartWrapper>
        {Object.entries(cellTypes).map(([tissue, cellTypeSummaries]) => {
          return (
            <Chart
              isScaled={isScaled}
              key={tissue}
              tissue={tissue}
              cellTypes={cellTypeSummaries}
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
