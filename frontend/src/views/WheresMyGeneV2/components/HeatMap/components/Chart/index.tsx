import { Tooltip } from "@czi-sds/components";
import { capitalize } from "lodash";
import cloneDeep from "lodash/cloneDeep";
import debounce from "lodash/debounce";
import throttle from "lodash/throttle";
import {
  Dispatch,
  memo,
  SetStateAction,
  useCallback,
  useContext,
  useEffect,
  useMemo,
  useRef,
  useState,
} from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { EMPTY_ARRAY, noop } from "src/common/constants/utils";
import {
  CellTypeRow,
  COMPARE_OPTION_ID_FOR_AGGREGATED,
  getOptionIdFromCellTypeViewId,
} from "src/common/queries/wheresMyGene";
import { getCompareOptionNameById } from "src/views/WheresMyGene/common/constants";
import { StateContext } from "src/views/WheresMyGene/common/store";
import {
  CellType,
  CellTypeGeneExpressionSummaryData,
  CellTypeSummary,
  ChartProps,
  GeneExpressionSummary,
  Tissue,
  ViewId,
} from "src/views/WheresMyGene/common/types";

import {
  getAllSerializedCellTypeMetadata,
  getGeneNames,
  getHeatmapHeight,
  getHeatmapWidth,
  hyphenize,
} from "../../../../../WheresMyGene/components/HeatMap/utils";
import { StyledHeatmapChart, StyledTooltipTable, tooltipCss } from "./style";

import {
  dataToChartFormat,
  grid,
  symbolSize,
  useChartItemStyle,
} from "./utils";

const CHART_DATA_TO_AXIS_ENCODING = {
  x: "geneIndex",
  y: "cellTypeIndex",
};

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
  maxExpression: number;
}

const BASE_DEBOUNCE_MS = 200;

const MAX_DEBOUNCE_MS = 2 * 1000;

const TOOLTIP_THROTTLE_MS = 100;

let handleDotHoverAnalytic: NodeJS.Timeout;

export default memo(function Chart({
  cellTypes: dataRows,
  selectedGeneData = EMPTY_ARRAY,
  setIsLoading,
  tissue,
  scaledMeanExpressionMax,
  scaledMeanExpressionMin,
  isScaled,
  echartsRendererMode,
  setAllChartProps,
  chartProps,
  maxExpression,
}: Props): JSX.Element {
  const [currentIndices, setCurrentIndices] = useState([-1, -1]);
  const [cursorOffset, setCursorOffset] = useState([-1, -1]);

  const ref = useRef<HTMLDivElement | null>(null);

  const [heatmapWidth, setHeatmapWidth] = useState(
    getHeatmapWidth(selectedGeneData)
  );

  const [heatmapHeight, setHeatmapHeight] = useState(
    getHeatmapHeight(dataRows)
  );

  useEffect(() => {
    setIsLoading((isLoading) => {
      if (isLoading[tissue]) return isLoading;

      return { ...isLoading, [tissue]: true };
    });
  }, [dataRows, selectedGeneData, setIsLoading, tissue]);

  const handleChartMouseMove = useMemo(() => {
    return throttle((params, chart) => {
      const { offsetX, offsetY, event } = params;
      const { pageX, pageY } = event as MouseEvent;

      if (!chart) return;

      const pointInGrid = chart.convertFromPixel("grid", [offsetX, offsetY]);

      setCursorOffset([pageX, pageY]);
      if (pointInGrid) setCurrentIndices(pointInGrid);
    }, TOOLTIP_THROTTLE_MS);
  }, []);

  // Update heatmap size
  useEffect(() => {
    setHeatmapWidth(getHeatmapWidth(selectedGeneData));
    setHeatmapHeight(getHeatmapHeight(dataRows));
  }, [dataRows, selectedGeneData]);

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
    debouncedIntegrateCellTypesAndGenes(dataRows, selectedGeneData);
  }, [selectedGeneData, dataRows, debouncedIntegrateCellTypesAndGenes]);

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
        const result: ChartProps = {
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
          const newAllChartProps = { ...allChartProps };
          newAllChartProps[tissue] = result;

          return newAllChartProps;
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

  const handleMouseLeave = useCallback(() => {
    // Handles race condition when a timeout is set after clearing
    setTimeout(() => {
      clearTimeout(handleDotHoverAnalytic);
    }, 100);
  }, []);

  const [hoveredGeneIndex, hoveredDataRowIndex] = currentIndices;

  const { compare } = useContext(StateContext);

  const tooltipContent = useMemo(() => {
    clearTimeout(handleDotHoverAnalytic);

    if (!chartProps) return null;

    const { chartData } = chartProps;

    // A row can either be a cell type or a tissue
    const dataPoint = chartData.find(
      ({ geneIndex, cellTypeIndex: dataRowIndex }) =>
        geneIndex === hoveredGeneIndex && dataRowIndex === hoveredDataRowIndex
    );

    const dataRow = dataRows[hoveredDataRowIndex];
    const gene = selectedGeneData[hoveredGeneIndex];

    if (!dataPoint || !dataRow || !gene) return null;

    const optionId = getOptionIdFromCellTypeViewId(
      dataPoint.id.split("-")[0] as ViewId
    );

    const { expressedCellCount } = dataPoint;

    const percentage = Number(((dataPoint.percentage || 0) * 100).toFixed(2));

    const tissuePercentage = Number(
      ((dataPoint.tissuePercentage || 0) * 100).toFixed(2)
    );

    const totalCellCount = dataRow.total_count;

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

    // cellTypeName will be the UBERON ID if the cell type is a tissue
    const cellTypeName = dataRow.cellTypeName.startsWith("UBERON")
      ? dataRow.name
      : dataRow.cellTypeName;

    const secondPanel = {
      dataRows: [
        { label: "Cell Type", value: cellTypeName },
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
          value: capitalize(dataRow.name),
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
    dataRows,
    hoveredGeneIndex,
    hoveredDataRowIndex,
    selectedGeneData,
    compare,
  ]);

  const tooltipClasses = useMemo(() => ({ tooltip: tooltipCss }), []);

  const tooltipPopperProps = useMemo(() => {
    return {
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
    };
  }, [cursorOffset]);

  const handleMouseMove = useCallback(() => {
    clearInterval(handleDotHoverAnalytic);
    if (tooltipContent?.props.data) {
      handleDotHoverAnalytic = setTimeout(() => {
        track(EVENTS.WMG_HEATMAP_DOT_HOVER);
      }, 2 * 1000);
    }
  }, [tooltipContent]);

  const { chartData, cellTypeMetadata, geneNames } = chartProps || {};

  const itemStyle = useChartItemStyle(isScaled, maxExpression);

  return (
    <Tooltip
      width="wide"
      classes={tooltipClasses}
      title={tooltipContent || <>No data</>}
      leaveDelay={0}
      placement="right-end"
      onMouseMove={handleMouseMove}
      PopperProps={tooltipPopperProps}
    >
      <StyledHeatmapChart
        height={heatmapHeight}
        width={heatmapWidth}
        ref={ref}
        id={`${hyphenize(tissue)}-chart`}
        onMouseLeave={handleMouseLeave}
        data={chartData}
        yAxisData={cellTypeMetadata}
        xAxisData={geneNames}
        echartsRendererMode={echartsRendererMode}
        onChartMouseMove={handleChartMouseMove}
        encode={CHART_DATA_TO_AXIS_ENCODING}
        itemStyle={itemStyle}
        symbolSize={symbolSize}
        grid={grid}
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
