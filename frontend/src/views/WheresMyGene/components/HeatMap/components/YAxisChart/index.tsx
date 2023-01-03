import { init } from "echarts";
import { memo, useContext, useEffect, useMemo, useRef, useState } from "react";
import { EMPTY_OBJECT, noop } from "src/common/constants/utils";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import Image from "next/image";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import { resetTissueCellTypes } from "src/views/WheresMyGene/common/store/actions";
import { CellType, Tissue } from "src/views/WheresMyGene/common/types";
import { useDeleteGenesAndCellTypes } from "../../hooks/useDeleteGenesAndCellTypes";
import { useUpdateYAxisChart } from "../../hooks/useUpdateYAxisChart";
import {
  CellTypeMetadata,
  deserializeCellTypeMetadata,
  getAllSerializedCellTypeMetadata,
  getHeatmapHeight,
} from "../../utils";
import ReplaySVG from "./icons/replay.svg";
import InfoSVG from "./icons/info-sign-icon.svg";
import {
  Container,
  ResetImageWrapper,
  TissueName,
  TissueWrapper,
  Wrapper,
  StyledImage,
  InfoButtonWrapper,
  Y_AXIS_TISSUE_WIDTH_PX,
} from "./style";

// the horizontal offset of the y-axis labels
// determined by echarts
const Y_AXIS_SVG_X_OFFSET = 12;
interface Props {
  cellTypes?: CellType[];
  hasDeletedCellTypes: boolean;
  availableCellTypes: CellType[];
  tissue: Tissue;
  tissueID: string;
  generateMarkerGenes: (cellType: CellType, tissueID: string) => void;
  selectedOrganismId: string;
}

type Coord = [number, number];
export default memo(function YAxisChart({
  cellTypes = [],
  hasDeletedCellTypes,
  availableCellTypes,
  tissue,
  tissueID,
  generateMarkerGenes,
  selectedOrganismId,
}: Props): JSX.Element {
  const tissueKey = tissue.replace(/\s+/g, "-");
  const yAxisSelectorQuery = `[data-test-id=cell-type-labels-${tissueKey}]`;

  const dispatch = useContext(DispatchContext);
  const { selectedTissues } = useContext(StateContext);

  const [isChartInitialized, setIsChartInitialized] = useState(false);
  const [yAxisInfoCoords, setYAxisInfoCoords] = useState<Coord[] | null>(null);
  const { cellTypeIdsToDelete, handleCellTypeClick } =
    useDeleteGenesAndCellTypes();

  const [yAxisChart, setYAxisChart] = useState<echarts.ECharts | null>(null);
  const yAxisRef = useRef(null);

  const [heatmapHeight, setHeatmapHeight] = useState(
    getHeatmapHeight(cellTypes)
  );
  const isMarkerGenes = get(FEATURES.MARKER_GENES) === BOOLEAN.TRUE;

  const setInfoCoordinates = () => {
    const topTissueKey = selectedTissues[0].replace(/\s+/g, "-");
    const containerTop = document.querySelector(
      `[data-test-id=cell-type-labels-${topTissueKey}]`
    );
    const container = document.querySelector(yAxisSelectorQuery);
    const textElements = container?.querySelectorAll(
      `text[transform*="translate(${Y_AXIS_SVG_X_OFFSET}"]`
    );
    if (container && containerTop && textElements) {
      const { left, top } = containerTop.getBoundingClientRect();
      const xOffset = left - Y_AXIS_TISSUE_WIDTH_PX;
      const yOffset = top;
      const infoCoords: Coord[] = [];
      textElements.forEach((el) => {
        const { right, top } = el.getBoundingClientRect();
        infoCoords.push([right - xOffset, top - yOffset]);
      });
      setYAxisInfoCoords(infoCoords);
    }
  };
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
    generateMarkerGenes,
    tissueID,
    availableCellTypes,
    yAxisChart,
    selectedOrganismId,
    isMarkerGenes,
    yAxisRef,
    yAxisInfoCoords,
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

  useEffect(() => {
    const targetNode = document.querySelector(yAxisSelectorQuery);
    const config = { attributes: true, childList: true, subtree: true };
    const callback = (mutationList: MutationRecord[]) => {
      for (const mutation of mutationList) {
        if (mutation.type === "childList" && mutation.target.nodeName === "g") {
          setInfoCoordinates();
          break;
        }
      }
    };
    const observer = new MutationObserver(callback);
    if (targetNode) {
      observer.observe(targetNode, config);
    }

    return () => {
      observer.disconnect();
    };
  }, [selectedTissues]);
  return (
    <Wrapper id={`${tissue}-y-axis`}>
      <TissueWrapper height={heatmapHeight}>
        <TissueName>{capitalize(tissue)}</TissueName>
        {hasDeletedCellTypes && (
          <ResetImageWrapper
            data-test-id="reset-cell-types"
            onClick={() => handleResetTissue(tissue)}
          >
            <Image
              src={ReplaySVG.src}
              width="12"
              height="12"
              alt="reset tissue cell types"
            />
          </ResetImageWrapper>
        )}
      </TissueWrapper>
      <Container
        data-test-id={`cell-type-labels-${tissueKey}`}
        height={heatmapHeight}
        ref={yAxisRef}
      />
      {yAxisInfoCoords &&
        isMarkerGenes &&
        yAxisInfoCoords.map((coord, i) => {
          const content = cellTypeMetadata[i];
          const cellType = isMarkerGenes ? deserializeCellTypeMetadata(content) : null;
          return (
            <InfoButtonWrapper
              id={content}
              key={content}
              left={coord[0]}
              top={coord[1]}
              onClick={() => {
                if (cellType) {
                  generateMarkerGenes(cellType, tissueID);
                }
              }}
              style={cellType && cellType.total_count < 25 ? {visibility: "hidden"} : {}}
            >
              <StyledImage
                id={"marker-gene-button"}
                src={InfoSVG.src}
                width="10"
                height="10"
                alt={`display marker genes for ${
                  deserializeCellTypeMetadata(content).name
                }`}
              />
            </InfoButtonWrapper>
          );
        })}
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
