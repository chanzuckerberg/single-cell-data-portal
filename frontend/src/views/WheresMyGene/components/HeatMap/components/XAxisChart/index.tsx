import { init } from "echarts";
import { useEffect, useRef, useState } from "react";
import { EMPTY_OBJECT } from "src/common/constants/utils";
import { useDeleteGenesAndCellTypes } from "../../hooks/useDeleteGenesAndCellTypes";
import { useUpdateXAxisChart } from "../../hooks/useUpdateXAxisChart";
import { getHeatmapWidth } from "../../utils";
import { XAxisContainer, XAxisMask, XAxisWrapper } from "./style";

interface Props {
  geneNames: string[];
}

export default function XAxisChart({ geneNames }: Props): JSX.Element {
  const [isChartInitialized, setIsChartInitialized] = useState(false);
  const [xAxisChart, setXAxisChart] = useState<echarts.ECharts | null>(null);
  const [heatmapWidth, setHeatmapWidth] = useState(getHeatmapWidth(geneNames));

  const { genesToDelete, handleGeneClick } = useDeleteGenesAndCellTypes();

  const xAxisRef = useRef(null);

  // Update heatmap size
  useEffect(() => {
    setHeatmapWidth(getHeatmapWidth(geneNames));
  }, [geneNames]);

  // Initialize charts
  useEffect(() => {
    const { current: xAxisCurrent } = xAxisRef;

    if (!xAxisCurrent || isChartInitialized) {
      return;
    }

    setIsChartInitialized(true);

    const xAxisChart = init(xAxisCurrent, EMPTY_OBJECT, {
      renderer: "svg",
      useDirtyRect: true,
    });

    setXAxisChart(xAxisChart);
  }, [xAxisRef, isChartInitialized]);

  // Bind xAxisChart events
  useEffect(() => {
    xAxisChart?.on("click", function (params) {
      /**
       * `value` is set by utils.getGeneNames()
       */
      const { value } = params;
      handleGeneClick(value as string);
    });
  }, [handleGeneClick, xAxisChart]);

  useUpdateXAxisChart({
    geneNames,
    genesToDelete,
    heatmapWidth,
    xAxisChart,
  });

  return (
    <XAxisWrapper width={heatmapWidth}>
      {/* (thuang): The extra div is needed to implement the mask */}
      <div>
        <XAxisContainer width={heatmapWidth} ref={xAxisRef} />
        <XAxisMask />
      </div>
    </XAxisWrapper>
  );
}
