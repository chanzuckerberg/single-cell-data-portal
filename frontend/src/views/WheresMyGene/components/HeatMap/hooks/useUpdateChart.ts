import throttle from "lodash/throttle";
import { useEffect, useMemo } from "react";
import { createChartOptions } from "../utils";
import { UPDATE_THROTTLE_MS } from "./common/constants";
import { ChartProps } from "./common/types";

interface Props {
  chart: echarts.ECharts | null;
  chartProps: ChartProps | null;
  isScaled: boolean;
}

export function useUpdateChart({ chart, chartProps, isScaled }: Props): void {
  const throttledUpdateChart = useMemo(() => {
    return throttle(
      ({ chart, chartProps, isScaled }: Props) => {
        if (!chart || !chartProps) {
          return;
        }

        const { cellTypeMetadata, geneNames, chartData } = chartProps;

        // (thuang): resize() needs to be called before setOption() to prevent
        // TypeError: Cannot read properties of undefined (reading 'shouldBePainted')
        chart.resize();

        chart.setOption(
          createChartOptions({
            cellTypeMetadata,
            chartData,
            geneNames,
            isScaled,
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
    return () => throttledUpdateChart.cancel();
  }, [throttledUpdateChart]);

  // Update the charts
  useEffect(() => {
    throttledUpdateChart({
      chart,
      chartProps,
      isScaled,
    });
  }, [chart, chartProps, throttledUpdateChart, isScaled]);
}
