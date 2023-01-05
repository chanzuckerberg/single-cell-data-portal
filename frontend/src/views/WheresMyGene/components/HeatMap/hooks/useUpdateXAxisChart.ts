import throttle from "lodash/throttle";
import { useEffect, useMemo } from "react";
import { createXAxisOptions } from "../utils";
import { UPDATE_THROTTLE_MS } from "./common/constants";

interface Props {
  geneNames?: string[];
  genesToDelete: string[];
  heatmapWidth: number;
  xAxisChart: echarts.ECharts | null;
}

export function useUpdateXAxisChart({
  geneNames,
  genesToDelete,
  heatmapWidth,
  xAxisChart,
}: Props): void {
  const throttledUpdateXAxisChart = useMemo(() => {
    return throttle(
      ({ geneNames = [], genesToDelete, xAxisChart, heatmapWidth }: Props) => {
        if (!xAxisChart) {
          return;
        }

        // (thuang): resize() needs to be called before setOption() to prevent
        // TypeError: Cannot read properties of undefined (reading 'shouldBePainted')
        xAxisChart.resize();

        xAxisChart.setOption(
          createXAxisOptions({
            geneNames,
            genesToDelete,
            heatmapWidth,
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
    return () => throttledUpdateXAxisChart.cancel();
  }, [throttledUpdateXAxisChart]);

  // Update xAxis chart
  useEffect(() => {
    throttledUpdateXAxisChart({
      geneNames,
      genesToDelete,
      heatmapWidth,
      xAxisChart,
    });
  }, [
    geneNames,
    genesToDelete,
    // (thuang): `heatmapWidth` is needed to make sure the chart resizes AFTER    
    // the DOM width has been updated
    heatmapWidth,
    xAxisChart,
    throttledUpdateXAxisChart,
  ]);
}
