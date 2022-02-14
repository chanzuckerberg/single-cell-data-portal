import { EChartsOption } from "echarts";
import throttle from "lodash/throttle";
import { useEffect, useMemo } from "react";
import { createChartOptions } from "../utils";
import { UPDATE_THROTTLE_MS } from "./common/constants";
import { ChartProps } from "./common/types";

export function useUpdateChart({
  chart,
  chartProps,
  commonOptions,
  isEchartGLAvailable,
}: {
  chart: echarts.ECharts | null;
  chartProps: ChartProps | null;
  commonOptions: EChartsOption;
  isEchartGLAvailable: boolean;
}): void {
  const throttledUpdateChart = useMemo(() => {
    return throttle(
      ({ chart, chartProps, isEchartGLAvailable, commonOptions }) => {
        if (!chart || !chartProps || !isEchartGLAvailable) {
          return;
        }

        const { cellTypeNames, geneNames } = chartProps;

        // (thuang): resize() needs to be called before setOption() to prevent
        // TypeError: Cannot read properties of undefined (reading 'shouldBePainted')
        chart.resize();

        chart.setOption(
          createChartOptions({ cellTypeNames, commonOptions, geneNames })
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
      commonOptions,
      isEchartGLAvailable,
    });
  }, [
    chart,
    commonOptions,
    isEchartGLAvailable,
    chartProps,
    throttledUpdateChart,
  ]);
}
