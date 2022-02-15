import { Intent, Spinner } from "@blueprintjs/core";
import { connect, DatasetComponentOption, EChartsOption, init } from "echarts";
import debounce from "lodash/debounce";
import { useContext, useEffect, useMemo, useRef, useState } from "react";
import { EMPTY_OBJECT, noop } from "src/common/constants/utils";
import { DispatchContext, State } from "../../common/store";
import { resetTissueCellTypes } from "../../common/store/actions";
import { CellTypeSummary } from "../../common/types";
import { ChartProps } from "./hooks/common/types";
import { useDeleteGenesAndCellTypes } from "./hooks/useDeleteGenesAndCellTypes";
import { useUpdateChart } from "./hooks/useUpdateChart";
import { useUpdateXAxisChart } from "./hooks/useUpdateXAxisChart";
import { useUpdateYAxisChart } from "./hooks/useUpdateYAxisChart";
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
  CellTypeMetadata,
  checkIsTissue,
  dataToChartFormat,
  deserializeCellTypeMetadata,
  getAllSerializedCellTypeMetadata,
  getGeneNames,
  getHeatmapHeight,
  getHeatmapWidth,
} from "./utils";

const DEBOUNCE_MS = 2 * 1000;

interface Props {
  cellTypes: CellTypeSummary[];
  data: CellTypeSummary[];
  genes: State["selectedGenes"];
  tissuesWithDeletedCellTypes: string[];
  allTissueCellTypes: { [region: string]: CellTypeSummary[] };
}

const ELEMENT_ID = "heat-map";

let isChartInitialized = false;

export default function HeatMap({
  cellTypes,
  data,
  genes,
  tissuesWithDeletedCellTypes,
  allTissueCellTypes,
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

  const {
    cellTypeIdsToDelete,
    genesToDelete,
    handleCellTypeClick,
    handleGeneClick,
  } = useDeleteGenesAndCellTypes();

  const dispatch = useContext(DispatchContext);

  // Loading state
  const [isLoading, setIsLoading] = useState(false);

  useEffect(() => {
    setIsLoading(true);
  }, [genes, cellTypes]);

  const debouncedDataToChartFormat = useMemo(() => {
    return debounce(
      (
        data: CellTypeSummary[],
        cellTypes: CellTypeSummary[],
        genes: State["selectedGenes"]
      ) => {
        const result = {
          cellTypeNames: getAllSerializedCellTypeMetadata(cellTypes),
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

  // Cancel debounce when unmounting
  useEffect(() => {
    return () => debouncedDataToChartFormat.cancel();
  }, [debouncedDataToChartFormat]);

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

  const [, setHandleYAxisChartClick] = useState(
    () => noop as (params: { value: CellTypeMetadata }) => void
  );

  // Bind yAxisChart events
  useEffect(() => {
    setHandleYAxisChartClick(
      (oldHandle: (params: { value: CellTypeMetadata }) => void) => {
        yAxisChart?.off("click", oldHandle);

        yAxisChart?.on("click", newHandle as never);

        return newHandle;
      }
    );

    function newHandle(params: { value: CellTypeMetadata }) {
      /**
       * `value` is set by utils.getAllSerializedCellTypeMetadata()
       */
      const { value } = params;

      if (checkIsTissue(value)) {
        if (!dispatch) return;

        const { tissue } = deserializeCellTypeMetadata(
          value as CellTypeMetadata
        );
        dispatch(resetTissueCellTypes(tissue, allTissueCellTypes[tissue]));
      } else {
        handleCellTypeClick(value);
      }
    }
  }, [
    setHandleYAxisChartClick,
    handleCellTypeClick,
    dispatch,
    allTissueCellTypes,
    yAxisChart,
  ]);

  const commonOptions = useMemo(() => {
    if (!chartProps) return {};

    const { chartData } = chartProps;

    return {
      animation: false,
      dataset: {
        source: chartData as DatasetComponentOption["source"],
      },
      hoverLayerThreshold: 10,
      progressive: 1e6,
    } as EChartsOption;
  }, [chartProps]);

  useUpdateChart({ chart, chartProps, commonOptions, isEchartGLAvailable });

  useUpdateXAxisChart({
    chartProps,
    commonOptions,
    genesToDelete,
    isEchartGLAvailable,
    xAxisChart,
  });

  useUpdateYAxisChart({
    cellTypeIdsToDelete,
    chartProps,
    commonOptions,
    isEchartGLAvailable,
    tissuesWithDeletedCellTypes,
    yAxisChart,
  });

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
