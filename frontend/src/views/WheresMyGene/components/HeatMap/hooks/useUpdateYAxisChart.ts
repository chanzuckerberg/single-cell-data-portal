import { ECharts } from "echarts";
import throttle from "lodash/throttle";
import { useEffect, useMemo } from "react";
import { CellTypeMetadata, createYAxisOptions } from "../utils";
import { UPDATE_THROTTLE_MS } from "./common/constants";

interface Props {
  cellTypeIdsToDelete: string[];
  cellTypeMetadata: CellTypeMetadata[];
  heatmapHeight: number;
  yAxisChart: ECharts | null;
}

export function useUpdateYAxisChart({
  cellTypeIdsToDelete,
  cellTypeMetadata,
  heatmapHeight,
  yAxisChart,
}: Props): void {
  const throttledUpdateYAxisChart = useMemo(() => {
    return throttle(
      ({ cellTypeIdsToDelete, cellTypeMetadata = [], yAxisChart }: Props) => {
        if (!cellTypeMetadata.length || !yAxisChart) {
          return;
        }

        // (thuang): resize() needs to be called before setOption() to prevent
        // TypeError: Cannot read properties of undefined (reading 'shouldBePainted')
        yAxisChart.resize();

        yAxisChart.setOption(
          createYAxisOptions({
            cellTypeIdsToDelete,
            cellTypeMetadata,
          })
        );
      },
      UPDATE_THROTTLE_MS,
      // (thuang): Trailing guarantees that the last call to the function will
      // be executed
      { trailing: true }
    );
  }, []);

  useEffect(() => {
    return () => throttledUpdateYAxisChart.cancel();
  }, [throttledUpdateYAxisChart]);

  // Update the yAxis chart
  useEffect(() => {
    throttledUpdateYAxisChart({
      cellTypeIdsToDelete,
      cellTypeMetadata,
      heatmapHeight,
      yAxisChart,
    });
  }, [
    cellTypeIdsToDelete,
    cellTypeMetadata,
    // (thuang): `heatmapHeight` is needed to make sure the chart resizes AFTER
    // the DOM height has been updated
    heatmapHeight,
    yAxisChart,
    throttledUpdateYAxisChart,
  ]);
}
