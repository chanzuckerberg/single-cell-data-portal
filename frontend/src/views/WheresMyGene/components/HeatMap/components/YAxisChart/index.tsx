import { init } from "echarts";
import Image from "next/image";
import { memo, useContext, useEffect, useMemo, useRef, useState } from "react";
import { EMPTY_OBJECT, noop } from "src/common/constants/utils";
import { DispatchContext } from "src/views/WheresMyGene/common/store";
import { resetTissueCellTypes } from "src/views/WheresMyGene/common/store/actions";
import { CellType, Tissue } from "src/views/WheresMyGene/common/types";
import { useDeleteGenesAndCellTypes } from "../../hooks/useDeleteGenesAndCellTypes";
import { useUpdateYAxisChart } from "../../hooks/useUpdateYAxisChart";
import {
  CellTypeMetadata,
  getAllSerializedCellTypeMetadata,
  getHeatmapHeight,
} from "../../utils";
import ReplaySVG from "./icons/replay.svg";
import {
  Container,
  ResetImageWrapper,
  TissueName,
  TissueWrapper,
  Wrapper,
} from "./style";

interface Props {
  cellTypes?: CellType[];
  hasDeletedCellTypes: boolean;
  availableCellTypes: CellType[];
  tissue: Tissue;
}

export default memo(function YAxisChart({
  cellTypes = [],
  hasDeletedCellTypes,
  availableCellTypes,
  tissue,
}: Props): JSX.Element {
  const dispatch = useContext(DispatchContext);

  const [isChartInitialized, setIsChartInitialized] = useState(false);

  const { cellTypeIdsToDelete, handleCellTypeClick } =
    useDeleteGenesAndCellTypes();

  const [yAxisChart, setYAxisChart] = useState<echarts.ECharts | null>(null);
  const yAxisRef = useRef(null);

  const [heatmapHeight, setHeatmapHeight] = useState(
    getHeatmapHeight(cellTypes)
  );

  // Initialize charts
  useEffect(() => {
    const { current: yAxisCurrent } = yAxisRef;

    if (!yAxisCurrent || isChartInitialized) {
      return;
    }

    setIsChartInitialized(true);

    const yAxisChart = init(yAxisCurrent, EMPTY_OBJECT, {
      renderer: "svg",
      useDirtyRect: true,
    });

    setYAxisChart(yAxisChart);
  }, [isChartInitialized]);

  // Update heatmap size
  useEffect(() => {
    setHeatmapHeight(getHeatmapHeight(cellTypes));
  }, [cellTypes]);

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

      handleCellTypeClick(value);
    }
  }, [
    setHandleYAxisChartClick,
    handleCellTypeClick,
    dispatch,
    availableCellTypes,
    yAxisChart,
  ]);

  const cellTypeMetadata = useMemo(() => {
    return getAllSerializedCellTypeMetadata(cellTypes, tissue);
  }, [cellTypes, tissue]);

  useUpdateYAxisChart({
    cellTypeIdsToDelete,
    cellTypeMetadata,
    heatmapHeight,
    yAxisChart,
  });

  return (
    <Wrapper>
      <TissueWrapper height={heatmapHeight}>
        <TissueName>{capitalize(tissue)}</TissueName>
        {hasDeletedCellTypes && (
          <ResetImageWrapper onClick={() => handleResetTissue(tissue)}>
            <Image
              src={ReplaySVG.src}
              width="12"
              height="12"
              alt="reset tissue cell types"
            />
          </ResetImageWrapper>
        )}
      </TissueWrapper>
      <Container height={heatmapHeight} ref={yAxisRef} />
    </Wrapper>
  );

  function handleResetTissue(tissue: Tissue) {
    if (!dispatch) return;

    dispatch(resetTissueCellTypes(tissue, availableCellTypes));
  }
});

function capitalize(str: string): string {
  return str.charAt(0).toUpperCase() + str.slice(1);
}
