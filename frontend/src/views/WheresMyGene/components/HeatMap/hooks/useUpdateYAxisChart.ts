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
        if (!yAxisChart) {
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
        // (alec): this is a hack to get echarts to display text overflow.
        const el = document.getElementById("cell-type-labels-axis");

        // @ts-ignore: style definitely exists on element but maybe it is not typed?
        el.children[0].style.overflow = "visible";
        // @ts-ignore: style definitely exists on element but maybe it is not typed?
        el.children[0].children[0].style.overflow = "visible";
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
