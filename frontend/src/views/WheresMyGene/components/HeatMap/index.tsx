import { Intent, Spinner } from "@blueprintjs/core";
import { connect, DatasetComponentOption, EChartsOption, init } from "echarts";
import debounce from "lodash/debounce";
import { useEffect, useMemo, useRef, useState } from "react";
import { EMPTY_OBJECT } from "src/common/constants/utils";
import { CellTypeAndGenes, Gene } from "../../common/types";
import {
  ChartContainer,
  Container,
  Loader,
  XAxisContainer,
  XAxisMask,
  XAxisWrapper,
  YAxisContainer,
} from "./style";
import {
  ChartFormat,
  createChartOptions,
  createXAxisOptions,
  createYAxisOptions,
  dataToChartFormat,
  getCellTypeNames,
  getGeneNames,
  getHeatmapHeight,
  getHeatmapWidth,
} from "./utils";

interface ChartProps {
  chartData: ChartFormat[];
  geneNames: string[];
  cellTypeNames: string[];
}

const DEBOUNCE_MS = 2 * 1000;

interface Props {
  cellTypes: CellTypeAndGenes[];
  data: CellTypeAndGenes[];
  genes: Gene[];
}

const ELEMENT_ID = "heat-map";

let isChartInitialized = false;

export default function HeatMap({
  cellTypes,
  data,
  genes,
}: Props): JSX.Element {
  const [chart, setChart] = useState<echarts.ECharts | null>(null);
  const [xAxisChart, setXAxisChart] = useState<echarts.ECharts | null>(null);
  const [yAxisChart, setYAxisChart] = useState<echarts.ECharts | null>(null);
  const [isEchartGLAvailable, setIsEchartGLAvailable] = useState(false);
  const [heatmapWidth, setHeatmapWidth] = useState(getHeatmapWidth(genes));
  const [heatmapHeight, setHeatmapHeight] = useState(
    getHeatmapHeight(cellTypes)
  );

  const [chartProps, setChartProps] = useState<ChartProps | null>(null);

  // Loading state
  const [isLoading, setIsLoading] = useState(false);

  useEffect(() => {
    setIsLoading(true);
  }, [genes]);

  const debouncedDataToChartFormat = useMemo(() => {
    return debounce(
      (data, cellTypes, genes) => {
        const result = {
          cellTypeNames: getCellTypeNames(cellTypes),
          chartData: dataToChartFormat(data, cellTypes, genes),
          geneNames: getGeneNames(genes),
        };

        setChartProps(result);

        setIsLoading(false);
      },
      DEBOUNCE_MS,
      { leading: false }
    );
  }, []);

  useEffect(() => {
    debouncedDataToChartFormat(data, cellTypes, genes);
  }, [data, cellTypes, genes, debouncedDataToChartFormat]);

  const ref = useRef(null);
  const xAxisRef = useRef(null);
  const yAxisRef = useRef(null);

  // Update heatmap size
  useEffect(() => {
    setHeatmapWidth(getHeatmapWidth(genes));
    setHeatmapHeight(getHeatmapHeight(cellTypes));
  }, [cellTypes, genes]);

  // Import echarts-gl
  useEffect(() => {
    importEchartsGl();

    async function importEchartsGl(): Promise<void> {
      await import("echarts-gl");
      setIsEchartGLAvailable(true);
    }
  }, []);

  // Initialize charts
  useEffect(() => {
    const { current } = ref;
    const { current: xAxisCurrent } = xAxisRef;
    const { current: yAxisCurrent } = yAxisRef;

    if (
      !current ||
      !xAxisCurrent ||
      !yAxisCurrent ||
      isChartInitialized ||
      !isEchartGLAvailable
    ) {
      return;
    }

    isChartInitialized = true;

    const chart = init(current, EMPTY_OBJECT, { useDirtyRect: true });
    const xAxisChart = init(xAxisCurrent, EMPTY_OBJECT, {
      useDirtyRect: true,
    });
    const yAxisChart = init(yAxisCurrent, EMPTY_OBJECT, {
      useDirtyRect: true,
    });

    connect([chart, xAxisChart, yAxisChart]);

    setChart(chart);
    setXAxisChart(xAxisChart);
    setYAxisChart(yAxisChart);
  }, [ref, xAxisRef, isEchartGLAvailable]);

  // Update the charts
  useEffect(() => {
    if (
      !chart ||
      !chartProps ||
      !xAxisChart ||
      !yAxisChart ||
      !isEchartGLAvailable
    ) {
      return;
    }

    const { chartData, cellTypeNames, geneNames } = chartProps;

    const commonOptions: EChartsOption = {
      animation: false,
      dataset: {
        source: chartData as DatasetComponentOption["source"],
      },
      hoverLayerThreshold: 10,
      progressive: 1e6,
    };

    chart.setOption(
      createChartOptions({ cellTypeNames, commonOptions, geneNames })
    );

    xAxisChart.setOption(
      createXAxisOptions({ cellTypeNames, commonOptions, geneNames })
    );

    yAxisChart.setOption(
      createYAxisOptions({ cellTypeNames, commonOptions, geneNames })
    );

    chart.resize();
    xAxisChart.resize();
    yAxisChart.resize();
  }, [chart, xAxisChart, yAxisChart, isEchartGLAvailable, chartProps]);

  return (
    <Container>
      {isLoading ? (
        <Loader>
          <Spinner intent={Intent.PRIMARY} size={20} />
          Loading...
        </Loader>
      ) : null}

      <XAxisWrapper width={heatmapWidth}>
        {/* (thuang): The extra div is needed to implement the mask */}
        <div>
          <XAxisContainer width={heatmapWidth} ref={xAxisRef} />
          <XAxisMask />
        </div>
      </XAxisWrapper>
      <YAxisContainer height={heatmapHeight} ref={yAxisRef} />
      <ChartContainer
        height={heatmapHeight}
        width={heatmapWidth}
        ref={ref}
        id={ELEMENT_ID}
      />
    </Container>
  );
}
