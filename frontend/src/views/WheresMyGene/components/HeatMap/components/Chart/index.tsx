import { Tooltip } from "czifui";
import { init } from "echarts";
import cloneDeep from "lodash/cloneDeep";
import throttle from "lodash/throttle";
import {
  Dispatch,
  memo,
  SetStateAction,
  useEffect,
  useMemo,
  useRef,
  useState,
} from "react";
import { EMPTY_ARRAY, EMPTY_OBJECT, noop } from "src/common/constants/utils";
import {
  CellType,
  CellTypeGeneExpressionSummaryData,
  CellTypeSummary,
  GeneExpressionSummary,
  Tissue,
} from "src/views/WheresMyGene/common/types";
import { ChartProps } from "../../hooks/common/types";
import { useUpdateChart } from "../../hooks/useUpdateChart";
import {
  dataToChartFormat,
  getAllSerializedCellTypeMetadata,
  getGeneNames,
  getHeatmapHeight,
  HEAT_MAP_BASE_CELL_WIDTH_PX
} from "../../utils";
import {
  ChartContainer,
  StyledTooltipTable,
  tooltipCss,
  Wrapper,
} from "./style";

interface Props {
  cellTypes: CellType[];
  setIsLoading: Dispatch<
    SetStateAction<{
      [tissue: Tissue]: boolean;
    }>
  >;
  tissue: Tissue;
  gene: GeneExpressionSummary;
  scaledMeanExpressionMax: number;
  scaledMeanExpressionMin: number;
  isScaled: boolean;
}

const TOOLTIP_THROTTLE_MS = 100;

export default memo(function Chart({
  cellTypes,
  setIsLoading,
  tissue,
  gene,
  scaledMeanExpressionMax,
  scaledMeanExpressionMin,
  isScaled,
}: Props): JSX.Element {
  const [currentIndices, setCurrentIndices] = useState([-1, -1]);
  const [cursorOffset, setCursorOffset] = useState([-1, -1]);

  const [isChartInitialized, setIsChartInitialized] = useState(false);

  const [chart, setChart] = useState<echarts.ECharts | null>(null);
  const ref = useRef(null);

  const heatmapWidth = HEAT_MAP_BASE_CELL_WIDTH_PX;
  
  const [heatmapHeight, setHeatmapHeight] = useState(
    getHeatmapHeight(cellTypes)
  );

  const [chartProps, setChartProps] = useState<ChartProps | null>(null);

  useEffect(() => {
    setIsLoading((isLoading) => ({ ...isLoading, [tissue]: true }));
  }, [cellTypes, gene, setIsLoading, tissue]);

  const throttledSetCurrentIndices = useMemo(() => {
    return throttle((params, chart) => {
      const { offsetX, offsetY, event } = params;
      const { pageX, pageY } = event;

      if (!chart) return;

      const pointInGrid = chart.convertFromPixel("grid", [offsetX, offsetY]);

      setCursorOffset([pageX, pageY]);
      setCurrentIndices(pointInGrid);
    }, TOOLTIP_THROTTLE_MS);
  }, []);

  // Initialize charts
  useEffect(() => {
    const { current } = ref;

    if (!current || isChartInitialized) {
      return;
    }

    setIsChartInitialized(true);

    const chart = init(current, EMPTY_OBJECT, { useDirtyRect: true });
    chart.getZr().on("mousemove", function (params) {
      throttledSetCurrentIndices(params, chart);
    });

    setChart(chart);
  }, [ref, isChartInitialized, throttledSetCurrentIndices]);

  // Update heatmap size
  useEffect(() => {
    setHeatmapHeight(getHeatmapHeight(cellTypes));
  }, [cellTypes]);

  useUpdateChart({ chart, chartProps, isScaled });

  // Calculate cellTypeSummaries
  /**
   * This is the formatted data that we use to render the heatmap.
   */
  const [cellTypeSummaries, setCellTypeSummaries] =
    useState<CellTypeSummary[]>(EMPTY_ARRAY);

  /**
   * Performance optimization:
   * We only format and `setCellTypeSummaries()` after the watch list has stopped changing for
   * `getDebounceMs()`
   */
  useEffect(() => {
    setCellTypeSummaries(
      integrateCelTypesAndGenes({
        cellTypes,
        geneExpressionSummaries: [gene],
      })
    );
  }, [gene, cellTypes]);

  // Generate chartProps

  useEffect(() => {
    setChartProps({
      cellTypeMetadata: getAllSerializedCellTypeMetadata(
        cellTypeSummaries,
        tissue
      ),
      chartData: dataToChartFormat({
        cellTypeSummaries,
        genes: [gene],
        scaledMeanExpressionMax,
        scaledMeanExpressionMin,
      }),
      geneNames: getGeneNames([gene]),
    });
    setIsLoading((isLoading) => ({ ...isLoading, [tissue]: false }));
  }, [cellTypeSummaries, gene, tissue, scaledMeanExpressionMax, scaledMeanExpressionMin, setIsLoading]);

  const [_, hoveredCellTypeIndex] = currentIndices;
  
  const tooltipContent = useMemo(() => {
    if (!chartProps) return null;

    const { chartData } = chartProps;
    const dataPoint = chartData.find(
      ({ id, cellTypeIndex }) =>
      {
        return id.split('-').at(-1) === gene.name && cellTypeIndex === hoveredCellTypeIndex
      }
    );
    const cellType = cellTypes[hoveredCellTypeIndex];

    if (!dataPoint || !cellType || !gene) return null;

    const { expressedCellCount } = dataPoint;

    const percentage = Number(((dataPoint.percentage || 0) * 100).toFixed(2));

    const tissuePercentage = Number(
      ((dataPoint.tissuePercentage || 0) * 100).toFixed(2)
    );

    const totalCellCount = Math.round((expressedCellCount / percentage) * 100);

    const data = [
      {
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
      },
      {
        dataRows: [
          { label: "Cell Type", value: cellType.name },
          {
            label: "Tissue Composition",
            value: tissuePercentage + "%" || "",
          },
        ],
      },
      {
        dataRows: [{ label: "Gene Symbol", value: gene.name || "" }],
      },
    ];

    return <StyledTooltipTable data={data || undefined} />;
  }, [
    chartProps,
    cellTypes,
    hoveredCellTypeIndex,
    gene
  ]);

  const tooltipClasses = useMemo(() => ({ tooltip: tooltipCss }), []);

  return (
    <Wrapper height={heatmapHeight} width={heatmapWidth}>
      <Tooltip
        width="wide"
        classes={tooltipClasses}
        title={tooltipContent || <>No data</>}
        leaveDelay={0}
        placement="right-end"
        PopperProps={{
          anchorEl: {
            clientHeight: 0,
            clientWidth: 0,
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
          modifiers: {
            offset: { offset: "0,20" },
          },
        }}
      >
        <ChartContainer height={heatmapHeight} width={heatmapWidth} ref={ref} />
      </Tooltip>
    </Wrapper>
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
  geneExpressionSummaries: GeneExpressionSummary[] | undefined;
}): CellTypeSummary[] {
  const geneMaps = geneExpressionSummaries.map((geneExpressionSummary) =>
    rawGeneDataToMap(geneExpressionSummary)
  );

  const cellTypeSummaries: CellTypeSummary[] = cloneDeep(cellTypes);

  return cellTypeSummaries.map((cellTypeSummary) => {
    const { id } = cellTypeSummary;

    for (const [name, geneMap] of geneMaps) {
      const columnData = geneMap.get(id);

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
    new Map(cellTypeGeneExpressionSummaries?.map((row) => [row.id, row])),
  ];
}