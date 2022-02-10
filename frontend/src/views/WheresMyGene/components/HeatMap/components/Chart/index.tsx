import { init } from "echarts";
import cloneDeep from "lodash/cloneDeep";
import debounce from "lodash/debounce";
import {
  Dispatch,
  memo,
  SetStateAction,
  useEffect,
  useMemo,
  useRef,
  useState,
} from "react";
import { EMPTY_ARRAY, EMPTY_OBJECT } from "src/common/constants/utils";
import {
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
  getHeatmapWidth,
} from "../../utils";
import { ChartContainer, Wrapper } from "./style";

interface Props {
  cellTypes: CellTypeSummary[];
  selectedGeneData: GeneExpressionSummary[];
  setIsLoading: Dispatch<
    SetStateAction<{
      lung: boolean;
    }>
  >;
  tissue: Tissue;
}

const DEBOUNCE_MS = 2 * 1000;

export default memo(function Chart({
  cellTypes,
  selectedGeneData,
  setIsLoading,
  tissue,
}: Props): JSX.Element {
  const [isChartInitialized, setIsChartInitialized] = useState(false);

  const [chart, setChart] = useState<echarts.ECharts | null>(null);
  const ref = useRef(null);

  const [heatmapWidth, setHeatmapWidth] = useState(
    getHeatmapWidth(selectedGeneData)
  );
  const [heatmapHeight, setHeatmapHeight] = useState(
    getHeatmapHeight(cellTypes)
  );

  const [chartProps, setChartProps] = useState<ChartProps | null>(null);

  // Initialize charts
  useEffect(() => {
    const { current } = ref;

    if (!current || isChartInitialized) {
      return;
    }

    setIsChartInitialized(true);

    const chart = init(current, EMPTY_OBJECT, { useDirtyRect: true });

    setChart(chart);
  }, [ref, isChartInitialized]);

  // Update heatmap size
  useEffect(() => {
    setHeatmapWidth(getHeatmapWidth(selectedGeneData));
    setHeatmapHeight(getHeatmapHeight(cellTypes));
  }, [cellTypes, selectedGeneData]);

  useUpdateChart({ chart, chartProps });

  // Calculate cellTypeSummaries
  /**
   * This is the formatted data that we use to render the heatmap.
   */
  const [cellTypeSummaries, setCellTypeSummaries] =
    useState<CellTypeSummary[]>(EMPTY_ARRAY);

  const debouncedIntegrateCellTypesAndGenes = useMemo(() => {
    return debounce(
      (cellTypes: CellTypeSummary[], geneData) => {
        setCellTypeSummaries(integrateCelTypesAndGenes(cellTypes, geneData));
      },
      DEBOUNCE_MS,
      { leading: false }
    );
  }, []);

  // Cancel debounce when unmounting
  useEffect(() => {
    return () => debouncedIntegrateCellTypesAndGenes.cancel();
  }, [debouncedIntegrateCellTypesAndGenes]);

  /**
   * Performance optimization:
   * We only format and `setCellTypeSummaries()` after the watch list has stopped changing for
   * `DEBOUNCE_MS`
   */
  useEffect(() => {
    debouncedIntegrateCellTypesAndGenes(cellTypes, selectedGeneData);
  }, [selectedGeneData, cellTypes, debouncedIntegrateCellTypesAndGenes]);

  // Generate chartProps
  const debouncedDataToChartFormat = useMemo(() => {
    return debounce(
      (
        cellTypeSummaries: CellTypeSummary[],
        selectedGeneData: GeneExpressionSummary[]
      ) => {
        const result = {
          cellTypeMetadata: getAllSerializedCellTypeMetadata(cellTypeSummaries),
          chartData: dataToChartFormat(cellTypeSummaries, selectedGeneData),
          geneNames: getGeneNames(selectedGeneData),
        };

        setChartProps(result);

        setIsLoading((isLoading) => {
          return {
            ...isLoading,
            [tissue]: false,
          };
        });
      },
      DEBOUNCE_MS,
      { leading: false }
    );
  }, [setIsLoading, tissue]);

  useEffect(() => {
    debouncedDataToChartFormat(cellTypeSummaries, selectedGeneData);
  }, [cellTypeSummaries, selectedGeneData, debouncedDataToChartFormat]);

  // Cancel debounce when unmounting
  useEffect(() => {
    return () => debouncedDataToChartFormat.cancel();
  }, [debouncedDataToChartFormat]);

  return (
    <Wrapper height={heatmapHeight} width={heatmapWidth}>
      <ChartContainer height={heatmapHeight} width={heatmapWidth} ref={ref} />
    </Wrapper>
  );
});

/**
 * Adds gene expressions to the cell types.
 */
function integrateCelTypesAndGenes(
  cellTypeSummaries: CellTypeSummary[],
  geneExpressionSummaries: GeneExpressionSummary[]
): CellTypeSummary[] {
  const geneMaps = geneExpressionSummaries.map((geneExpressionSummary) =>
    rawGeneDataToMap(geneExpressionSummary)
  );

  const newCellTypeSummaries = cloneDeep(cellTypeSummaries);

  return newCellTypeSummaries.map((cellTypeSummary) => {
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
  gene: GeneExpressionSummary
): [string, Map<string, CellTypeGeneExpressionSummaryData>] {
  const { cellTypeGeneExpressionSummaries, name } = gene;

  return [
    name,
    new Map(cellTypeGeneExpressionSummaries?.map((row) => [row.id, row])),
  ];
}
