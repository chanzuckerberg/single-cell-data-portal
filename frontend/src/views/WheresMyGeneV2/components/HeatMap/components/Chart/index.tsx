import { Tooltip } from "czifui";
import { ECharts, init } from "echarts";
import { capitalize } from "lodash";
import cloneDeep from "lodash/cloneDeep";
import debounce from "lodash/debounce";
import throttle from "lodash/throttle";
import {
  Dispatch,
  memo,
  SetStateAction,
  useContext,
  useEffect,
  useMemo,
  useRef,
  useState,
} from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { EMPTY_ARRAY, EMPTY_OBJECT, noop } from "src/common/constants/utils";
import {
  CellTypeRow,
  COMPARE_OPTION_ID_FOR_AGGREGATED,
  getOptionIdFromCellTypeViewId,
} from "src/common/queries/wheresMyGeneV2";
import { getCompareOptionNameById } from "src/views/WheresMyGeneV2/common/constants";
import { StateContext } from "src/views/WheresMyGeneV2/common/store";
import {
  CellType,
  CellTypeGeneExpressionSummaryData,
  CellTypeSummary,
  GeneExpressionSummary,
  Tissue,
  ViewId,
} from "src/views/WheresMyGeneV2/common/types";
import { ChartProps } from "../../hooks/common/types";
import { useUpdateChart } from "../../hooks/useUpdateChart";
import {
  dataToChartFormat,
  getAllSerializedCellTypeMetadata,
  getGeneNames,
  getHeatmapHeight,
  getHeatmapWidth,
} from "../../utils";
import { ChartContainer, StyledTooltipTable, tooltipCss } from "./style";

interface Props {
  cellTypes: CellTypeRow[];
  selectedGeneData?: (GeneExpressionSummary | undefined)[];
  setIsLoading: Dispatch<
    SetStateAction<{
      [tissue: Tissue]: boolean;
    }>
  >;
  tissue: Tissue;
  scaledMeanExpressionMax: number;
  scaledMeanExpressionMin: number;
  isScaled: boolean;
  echartsRendererMode: "svg" | "canvas";
  setAllChartProps: Dispatch<
    SetStateAction<{
      [tissue: string]: ChartProps;
    }>
  >;
  chartProps: ChartProps;
}

const BASE_DEBOUNCE_MS = 200;

const MAX_DEBOUNCE_MS = 2 * 1000;

const TOOLTIP_THROTTLE_MS = 100;

let handleDotHoverAnalytic: NodeJS.Timeout;

export default memo(function Chart({
  cellTypes,
  selectedGeneData = EMPTY_ARRAY,
  setIsLoading,
  tissue,
  scaledMeanExpressionMax,
  scaledMeanExpressionMin,
  isScaled,
  echartsRendererMode,
  setAllChartProps,
  chartProps,
}: Props): JSX.Element {
  const [currentIndices, setCurrentIndices] = useState([-1, -1]);
  const [cursorOffset, setCursorOffset] = useState([-1, -1]);

  const [isChartInitialized, setIsChartInitialized] = useState(false);

  const [chart, setChart] = useState<ECharts | null>(null);
  const ref = useRef<HTMLDivElement | null>(null);

  const [heatmapWidth, setHeatmapWidth] = useState(
    getHeatmapWidth(selectedGeneData)
  );

  const [heatmapHeight, setHeatmapHeight] = useState(
    getHeatmapHeight(cellTypes)
  );

  useEffect(() => {
    setIsLoading((isLoading) => ({ ...isLoading, [tissue]: true }));
  }, [cellTypes, selectedGeneData, setIsLoading, tissue]);

  const throttledSetCurrentIndices = useMemo(() => {
    return throttle((params, chart: ECharts) => {
      const { offsetX, offsetY, event } = params;
      const { pageX, pageY } = event;

      if (!chart) return;

      const pointInGrid = chart.convertFromPixel("grid", [offsetX, offsetY]);

      setCursorOffset([pageX, pageY]);
      if (pointInGrid) setCurrentIndices(pointInGrid);
    }, TOOLTIP_THROTTLE_MS);
  }, []);

  // Initialize charts
  useEffect(() => {
    const { current } = ref;

    if (
      !current ||
      isChartInitialized ||
      // (thuang): echart's `init()` will throw error if the container has 0 width or height
      current?.getAttribute("height") === "0" ||
      current?.getAttribute("width") === "0"
    ) {
      return;
    }

    setIsChartInitialized(true);

    const chart = init(current, EMPTY_OBJECT, {
      renderer: echartsRendererMode,
      useDirtyRect: true,
    });

    chart.getZr().on("mousemove", function (params) {
      throttledSetCurrentIndices(params, chart);
    });

    setChart(chart);
  }, [
    ref,
    isChartInitialized,
    throttledSetCurrentIndices,
    echartsRendererMode,
  ]);

  // Update heatmap size
  useEffect(() => {
    setHeatmapWidth(getHeatmapWidth(selectedGeneData));
    setHeatmapHeight(getHeatmapHeight(cellTypes));
  }, [cellTypes, selectedGeneData]);

  useUpdateChart({
    chart,
    chartProps,
    heatmapHeight,
    heatmapWidth,
    isScaled,
  });

  // Calculate cellTypeSummaries
  /**
   * This is the formatted data that we use to render the heatmap.
   */
  const [cellTypeSummaries, setCellTypeSummaries] =
    useState<CellTypeSummary[]>(EMPTY_ARRAY);

  const debouncedIntegrateCellTypesAndGenes = useMemo(() => {
    return debounce(
      (
        cellTypes: CellTypeRow[],
        geneData: Props["selectedGeneData"] = EMPTY_ARRAY
      ) => {
        setCellTypeSummaries(
          integrateCelTypesAndGenes({
            cellTypes,
            geneExpressionSummaries: geneData,
          })
        );
      },
      getDebounceMs(selectedGeneData.length),
      { leading: false }
    );
  }, [selectedGeneData]);

  // Cancel debounce when unmounting
  useEffect(() => {
    return () => debouncedIntegrateCellTypesAndGenes.cancel();
  }, [debouncedIntegrateCellTypesAndGenes]);

  /**
   * Performance optimization:
   * We only format and `setCellTypeSummaries()` after the watch list has stopped changing for
   * `getDebounceMs()`
   */
  useEffect(() => {
    debouncedIntegrateCellTypesAndGenes(cellTypes, selectedGeneData);
  }, [selectedGeneData, cellTypes, debouncedIntegrateCellTypesAndGenes]);

  // Generate chartProps
  const debouncedDataToChartFormat = useMemo(() => {
    return debounce(
      ({
        cellTypeSummaries,
        selectedGeneData = EMPTY_ARRAY,
        setAllChartProps,
      }: {
        cellTypeSummaries: CellTypeSummary[];
        selectedGeneData: Props["selectedGeneData"];
        setAllChartProps: Dispatch<
          SetStateAction<{
            [tissue: string]: ChartProps;
          }>
        >;
      }) => {
        const result = {
          cellTypeMetadata: getAllSerializedCellTypeMetadata(
            cellTypeSummaries,
            tissue
          ),
          chartData: dataToChartFormat({
            cellTypeSummaries,
            genes: selectedGeneData,
            scaledMeanExpressionMax,
            scaledMeanExpressionMin,
          }),
          geneNames: getGeneNames(selectedGeneData),
        };

        setAllChartProps((allChartProps) => {
          allChartProps[tissue] = result;
          return allChartProps;
        });

        setIsLoading((isLoading) => {
          if (!isLoading[tissue]) return isLoading;
          return { ...isLoading, [tissue]: false };
        });
      },
      getDebounceMs(selectedGeneData.length),
      { leading: false }
    );
  }, [
    selectedGeneData,
    setIsLoading,
    tissue,
    scaledMeanExpressionMax,
    scaledMeanExpressionMin,
  ]);

  useEffect(() => {
    debouncedDataToChartFormat({
      cellTypeSummaries,
      selectedGeneData,
      setAllChartProps,
    });
  }, [
    cellTypeSummaries,
    selectedGeneData,
    debouncedDataToChartFormat,
    setAllChartProps,
  ]);

  // Cancel debounce when unmounting
  useEffect(() => {
    return () => debouncedDataToChartFormat.cancel();
  }, [debouncedDataToChartFormat]);

  const [hoveredGeneIndex, hoveredCellTypeIndex] = currentIndices;

  const { compare } = useContext(StateContext);

  const tooltipContent = useMemo(() => {
    clearTimeout(handleDotHoverAnalytic);

    if (!chartProps) return null;

    const { chartData } = chartProps;

    const dataPoint = chartData.find(
      ({ geneIndex, cellTypeIndex }) =>
        geneIndex === hoveredGeneIndex && cellTypeIndex === hoveredCellTypeIndex
    );

    const cellType = cellTypes[hoveredCellTypeIndex];
    const gene = selectedGeneData[hoveredGeneIndex];

    if (!dataPoint || !cellType || !gene) return null;

    const optionId = getOptionIdFromCellTypeViewId(
      dataPoint.id.split("-")[0] as ViewId
    );

    const { expressedCellCount } = dataPoint;

    const percentage = Number(((dataPoint.percentage || 0) * 100).toFixed(2));

    const tissuePercentage = Number(
      ((dataPoint.tissuePercentage || 0) * 100).toFixed(2)
    );

    const totalCellCount = cellType.total_count;

    const firstPanel = {
      dataRows: [
        {
          label: "Expressed in Cells",
          value: `${percentage}% (${expressedCellCount} of ${totalCellCount} cells)`,
        },
        {
          label: "Gene Expression",
          value: (dataPoint.meanExpression || 0).toFixed(2),
        },
        {
          label: "Gene Expression, Scaled",
          value: (dataPoint.scaledMeanExpression || 0).toFixed(2),
        },
      ],
    };

    const secondPanel = {
      dataRows: [
        { label: "Cell Type", value: cellType.cellTypeName },
        {
          label: "Tissue Composition",
          value: tissuePercentage + "%" || "",
        },
      ],
    };

    if (compare && optionId !== COMPARE_OPTION_ID_FOR_AGGREGATED) {
      secondPanel.dataRows = [
        secondPanel.dataRows[0],
        {
          label: getCompareOptionNameById(compare),
          value: capitalize(cellType.name),
        },
        ...secondPanel.dataRows.slice(1),
      ];
    }

    const thirdPanel = {
      dataRows: [{ label: "Gene Symbol", value: gene.name || "" }],
    };

    const data = [firstPanel, secondPanel, thirdPanel];

    return <StyledTooltipTable data={data || undefined} />;
  }, [
    chartProps,
    cellTypes,
    hoveredGeneIndex,
    hoveredCellTypeIndex,
    selectedGeneData,
    compare,
  ]);

  const tooltipClasses = useMemo(() => ({ tooltip: tooltipCss }), []);

  return (
    <Tooltip
      width="wide"
      classes={tooltipClasses}
      title={tooltipContent || <>No data</>}
      leaveDelay={0}
      placement="right-end"
      onMouseMove={() => {
        clearInterval(handleDotHoverAnalytic);
        if (tooltipContent?.props.data) {
          handleDotHoverAnalytic = setTimeout(() => {
            track(EVENTS.WMG_HEATMAP_DOT_HOVER);
          }, 2 * 1000);
        }
      }}
      PopperProps={{
        anchorEl: {
          getBoundingClientRect: () => ({
            bottom: cursorOffset[1],
            height: 0,
            left: cursorOffset[0],
            right: cursorOffset[0],
            toJSON: noop,
            top: cursorOffset[1],
            width: 0,
            x: cursorOffset[0],
            y: cursorOffset[1],
          }),
        },
        modifiers: [
          {
            name: "offset",
            options: {
              offset: [0, 5],
            },
          },
        ],
      }}
    >
      <ChartContainer
        height={heatmapHeight}
        width={heatmapWidth}
        ref={ref}
        id={`${tissue.replace(/\s+/g, "-")}-chart`}
        onMouseLeave={() => {
          // Handles race condition when a timeout is set after clearing
          setTimeout(() => {
            clearTimeout(handleDotHoverAnalytic);
          }, 100);
        }}
      />
    </Tooltip>
  );
});

/**
 * Adds gene expressions to the selected cell types.
 */
function integrateCelTypesAndGenes({
  cellTypes,
  geneExpressionSummaries = EMPTY_ARRAY,
}: {
  cellTypes: CellType[];
  geneExpressionSummaries: Props["selectedGeneData"];
}): CellTypeSummary[] {
  const geneMaps = geneExpressionSummaries.map((geneExpressionSummary) =>
    rawGeneDataToMap(geneExpressionSummary)
  );

  const cellTypeSummaries: CellTypeSummary[] = cloneDeep(cellTypes);

  return cellTypeSummaries.map((cellTypeSummary) => {
    const { viewId } = cellTypeSummary;

    for (const [name, geneMap] of geneMaps) {
      const columnData = geneMap.get(viewId);

      if (columnData !== undefined) {
        cellTypeSummary.geneExpressions = {
          ...(cellTypeSummary.geneExpressions || {}),
          [name]: columnData,
        };
      }
    }

    return cellTypeSummary;
  });
}

function rawGeneDataToMap(
  gene?: GeneExpressionSummary
): [string, Map<string, CellTypeGeneExpressionSummaryData>] {
  if (!gene) return ["", new Map()];

  const { cellTypeGeneExpressionSummaries, name } = gene;

  return [
    name,
    new Map(cellTypeGeneExpressionSummaries?.map((row) => [row.viewId, row])),
  ];
}

const BROWSER_PARALLEL_CALL_LIMIT = 10;

function getDebounceMs(geneCount: number): number {
  if (geneCount <= BROWSER_PARALLEL_CALL_LIMIT) return 0;
  if (geneCount >= 100) return MAX_DEBOUNCE_MS;

  return Math.floor(geneCount / BROWSER_PARALLEL_CALL_LIMIT) * BASE_DEBOUNCE_MS;
}
