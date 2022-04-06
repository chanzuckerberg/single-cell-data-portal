import { LoadingIndicator } from "czifui";
import { memo, useEffect, useMemo, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
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

enum FirstLoadState {
  Initial,
  Loading,
  Loaded,
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
  const [firstLoad, setFirstLoad] = useState(FirstLoadState.Initial);

  // (thuang): We only want to send `WMG_HEATMAP_LOADED` event the first time it loads
  useEffect(() => {
    if (firstLoad === FirstLoadState.Loaded) return;
    if (firstLoad === FirstLoadState.Initial && isAnyTissueLoading(isLoading)) {
      setFirstLoad(FirstLoadState.Loading);
      return;
    }
    if (
      firstLoad === FirstLoadState.Loading &&
      !isAnyTissueLoading(isLoading)
    ) {
      track(EVENTS.WMG_HEATMAP_LOADED);
      setFirstLoad(FirstLoadState.Loaded);
    }
  }, [firstLoad, isLoading]);

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
        <LoaderComponent />
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

  function LoaderComponent() {
    return (
      <Loader>
        <LoadingIndicator sdsStyle="tag" />
      </Loader>
    );
  }
});

function isAnyTissueLoading(isLoading: { [tissue: Tissue]: boolean }) {
  return Object.values(isLoading).some((isLoading) => isLoading);
}

function setInitialIsLoading(cellTypes: Props["cellTypes"]) {
  return Object.keys(cellTypes).reduce((isLoading, tissue) => {
    return { ...isLoading, [tissue]: false };
  }, {});
}
