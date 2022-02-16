import { EChartsOption } from "echarts";
import throttle from "lodash/throttle";
import { useEffect, useMemo } from "react";
import { createXAxisOptions } from "../utils";
import { UPDATE_THROTTLE_MS } from "./common/constants";
import { ChartProps } from "./common/types";

export function useUpdateXAxisChart({
  chartProps,
  commonOptions,
  genesToDelete,
  isEchartGLAvailable,
  xAxisChart,
}: {
  chartProps: ChartProps | null;
  commonOptions: EChartsOption;
  genesToDelete: string[];
  isEchartGLAvailable: boolean;
  xAxisChart: echarts.ECharts | null;
}): void {
  const throttledUpdateXAxisChart = useMemo(() => {
    return throttle(
      ({
        chartProps,
        genesToDelete,
        isEchartGLAvailable,
        xAxisChart,
        commonOptions,
      }) => {
        if (!chartProps || !xAxisChart || !isEchartGLAvailable) {
          return;
        }

        const { cellTypeNames, geneNames } = chartProps;

        // (thuang): resize() needs to be called before setOption() to prevent
        // TypeError: Cannot read properties of undefined (reading 'shouldBePainted')
        xAxisChart.resize();

        xAxisChart.setOption(
          createXAxisOptions({
            cellTypeNames,
            commonOptions,
            geneNames,
            genesToDelete,
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
      chartProps,
      commonOptions,
      genesToDelete,
      isEchartGLAvailable,
      xAxisChart,
    });
  }, [
    chartProps,
    genesToDelete,
    isEchartGLAvailable,
    xAxisChart,
    commonOptions,
    throttledUpdateXAxisChart,
  ]);
}
