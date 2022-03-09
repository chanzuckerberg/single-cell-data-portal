import { Intent, Spinner } from "@blueprintjs/core";
import { memo, useMemo, useState } from "react";
import { State } from "../../common/store";
import {
  CellTypeSummary,
  GeneExpressionSummary,
  Tissue,
} from "../../common/types";
import Chart from "./components/Chart";
import XAxisChart from "./components/XAxisChart";
import YAxisChart from "./components/YAxisChart";
import { ChartWrapper, Container, Loader, YAxisWrapper } from "./style";
import { getHeatmapHeight, X_AXIS_CHART_HEIGHT_PX } from "./utils";

interface Props {
  cellTypes: { [tissue: Tissue]: CellTypeSummary[] };
  genes: State["selectedGenes"];
  tissuesWithDeletedCellTypes: string[];
  allTissueCellTypes: { [tissue: Tissue]: CellTypeSummary[] };
  selectedGeneData: GeneExpressionSummary[];
}

export default memo(function HeatMap({
  cellTypes,
  genes,
  tissuesWithDeletedCellTypes,
  allTissueCellTypes,
  selectedGeneData,
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
      {isAnyTissueLoading(isLoading) ? (
        <Loader>
          <Spinner intent={Intent.PRIMARY} size={20} />
          Loading...
        </Loader>
      ) : null}

      <XAxisChart geneNames={genes} />

      <YAxisWrapper height={yAxisWrapperHeight}>
        {Object.entries(cellTypes).map(([tissue, cellTypeSummaries]) => {
          return (
            <YAxisChart
              key={tissue}
              tissue={tissue}
              cellTypes={cellTypeSummaries}
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
              key={tissue}
              tissue={tissue}
              cellTypes={cellTypeSummaries}
              selectedGeneData={selectedGeneData}
              setIsLoading={setIsLoading}
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
