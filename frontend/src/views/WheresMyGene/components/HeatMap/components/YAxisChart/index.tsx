import { AnchorButton, Icon } from "@blueprintjs/core";
import { init } from "echarts";
import Image from "next/image";
import { memo, useContext, useEffect, useMemo, useRef, useState } from "react";
import { EMPTY_OBJECT, noop } from "src/common/constants/utils";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import { FetchMarkerGeneParams } from "src/common/queries/wheresMyGene";
import { DispatchContext } from "src/views/WheresMyGene/common/store";
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
} from "./style";

interface Props {
  cellTypes?: CellType[];
  hasDeletedCellTypes: boolean;
  availableCellTypes: CellType[];
  tissue: Tissue;
  tissueID: string;
  generateMarkerGenes: (args: FetchMarkerGeneParams) => void;
  selectedOrganismId: string;
}

export default memo(function YAxisChart({
  cellTypes = [],
  hasDeletedCellTypes,
  availableCellTypes,
  tissue,
  tissueID,
  generateMarkerGenes,
  selectedOrganismId,
}: Props): JSX.Element {
  const dispatch = useContext(DispatchContext);

  const [isChartInitialized, setIsChartInitialized] = useState(false);
  const [yAxisTextElements, setYAxisTextElements] = useState<NodeListOf<Element> | null>(null);
  const [elementToCellMetadata, setElementToCellMetadata] = useState<Map<Element, CellTypeMetadata>>(new Map());
  const { cellTypeIdsToDelete, handleCellTypeClick } =
    useDeleteGenesAndCellTypes();

  const [yAxisChart, setYAxisChart] = useState<echarts.ECharts | null>(null);
  const yAxisRef = useRef(null);

  const [heatmapHeight, setHeatmapHeight] = useState(
    getHeatmapHeight(cellTypes)
  );
  const isMarkerGenes = get(FEATURES.MARKER_GENES) === BOOLEAN.TRUE;

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

      if (isMarkerGenes) {
        const { id } = deserializeCellTypeMetadata(value);

        generateMarkerGenes({
          cellTypeID: id,
          organismID: selectedOrganismId,
          tissueID,
        });
      }
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
    yAxisRef
  ]);


  const cellTypeMetadata = useMemo(() => {
    return getAllSerializedCellTypeMetadata(cellTypes, tissue);
  }, [cellTypes, tissue]);

  useEffect(()=>{
    setTimeout(() => {
    const textElements = document.querySelector(`[data-test-id=cell-type-labels-${tissue}]`)?.querySelectorAll(`text[transform*="translate(12"]`)
    if (textElements) {
      setYAxisTextElements(textElements)
      const elToCell = new Map();
      textElements.forEach((el, i) => {
        elToCell.set(el, cellTypeMetadata[i]);
      })
      setElementToCellMetadata(elToCell);
    }
    
    }, 200)
  }, [yAxisChart, cellTypeMetadata])  

  useUpdateYAxisChart({
    cellTypeIdsToDelete,
    cellTypeMetadata,
    heatmapHeight,
    yAxisChart,
  });

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
        data-test-id={`cell-type-labels-${tissue}`}
        height={heatmapHeight}
        ref={yAxisRef}
      />
      {yAxisTextElements && Array.from(yAxisTextElements).map((el) => {
        const val = el.getBoundingClientRect().top-238
        const val2 = el.getBoundingClientRect().right
        const content = elementToCellMetadata.get(el)
        return (
          <div
            key={`${content ?? val}`}
            style={{
              position: "absolute",
              left: val2-278,
              top: val,
              cursor: "pointer",
            }}   
            onClick={() => {
              console.log(`Display marker genes for ${content}`)
            }}       
          >
            <Image
              src={InfoSVG.src}
              width="9"
              height="9"
              alt="display marker genes"
            />            
          </div>
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
