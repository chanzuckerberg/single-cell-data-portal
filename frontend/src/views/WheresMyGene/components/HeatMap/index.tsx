import { memo, useEffect, useMemo, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { State } from "../../common/store";
import { CellType, GeneExpressionSummary, Tissue } from "../../common/types";
import Loader from "../Loader";
import Chart from "./components/Chart";
import XAxisChart from "./components/XAxisChart";
import YAxisChart from "./components/YAxisChart";
import { ChartWrapper, Container, YAxisWrapper } from "./style";
import {
  getHeatmapHeight,
  HEAT_MAP_BASE_HEIGHT_PX,
  X_AXIS_CHART_HEIGHT_PX,
} from "./utils";

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

enum FirstLoadState {
  Initial,
  Loading,
  Loaded,
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

    // (thuang): We can't use Object.keys(cellTypes) here, because when genes
    // are not selected, we don't have cell types data from the API to populate
    // `cellTypes`
    const tissueCount = selectedTissues.length;

    const baseTissueHeightWhenNoGenesSelected =
      HEAT_MAP_BASE_HEIGHT_PX * tissueCount;

    return (
      Math.max(yAxisChartHeight, baseTissueHeightWhenNoGenesSelected) +
      X_AXIS_CHART_HEIGHT_PX
    );
  }, [cellTypes, selectedTissues]);

  return (
    <Container>
      {isLoadingAPI || isAnyTissueLoading(isLoading) ? <Loader /> : null}

      <XAxisChart geneNames={genes} />

      <YAxisWrapper height={yAxisWrapperHeight}>
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
      <ChartWrapper>
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
