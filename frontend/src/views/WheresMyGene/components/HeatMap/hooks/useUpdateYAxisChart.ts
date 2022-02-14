import { ECharts, EChartsOption } from "echarts";
import throttle from "lodash/throttle";
import { useEffect, useMemo } from "react";
import { createYAxisOptions } from "../utils";
import { UPDATE_THROTTLE_MS } from "./common/constants";
import { ChartProps } from "./common/types";

export function useUpdateYAxisChart({
  cellTypeIdsToDelete,
  chartProps,
  commonOptions,
  isEchartGLAvailable,
  tissuesWithDeletedCellTypes,
  yAxisChart,
}: {
  cellTypeIdsToDelete: string[];
  chartProps: ChartProps | null;
  commonOptions: EChartsOption;
  isEchartGLAvailable: boolean;
  tissuesWithDeletedCellTypes: string[];
  yAxisChart: ECharts | null;
}): void {
  const throttledUpdateYAxisChart = useMemo(() => {
    return throttle(
      ({
        cellTypeIdsToDelete,
        chartProps,
        isEchartGLAvailable,
        tissuesWithDeletedCellTypes,
        yAxisChart,
        commonOptions,
      }) => {
        if (!chartProps || !yAxisChart || !isEchartGLAvailable) {
          return;
        }

        const { cellTypeNames, geneNames } = chartProps;

        // (thuang): resize() needs to be called before setOption() to prevent
        // TypeError: Cannot read properties of undefined (reading 'shouldBePainted')
        yAxisChart.resize();

        yAxisChart.setOption(
          createYAxisOptions({
            cellTypeIdsToDelete,
            cellTypeNames,
            commonOptions,
            geneNames,
            tissuesWithDeletedCellTypes,
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
      chartProps,
      commonOptions,
      isEchartGLAvailable,
      tissuesWithDeletedCellTypes,
      yAxisChart,
    });
  }, [
    cellTypeIdsToDelete,
    chartProps,
    commonOptions,
    isEchartGLAvailable,
    tissuesWithDeletedCellTypes,
    yAxisChart,
    throttledUpdateYAxisChart,
  ]);
}
