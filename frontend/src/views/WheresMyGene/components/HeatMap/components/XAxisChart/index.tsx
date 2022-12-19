import { init } from "echarts";
import { useContext, useEffect, useRef, useState } from "react";
import { EMPTY_OBJECT } from "src/common/constants/utils";
import { removeSelectedGenes } from  "src/views/WheresMyGene/common/store/actions";
import { useDeleteGenesAndCellTypes } from "../../hooks/useDeleteGenesAndCellTypes";
import { useUpdateXAxisChart } from "../../hooks/useUpdateXAxisChart";
import { getHeatmapWidth } from "../../utils";
import { XAxisWrapper, XAxisContainer, GeneGroupWrapper, GeneGroupName, MarkerGeneHeader, MarkerGeneHeaderButton } from "./style";
import xSvg from "./icons/x-button.svg";
import Image from "next/image";
import { DispatchContext } from "src/views/WheresMyGene/common/store";


interface Props {
  geneNames: string[];
  leftOffset: number;
  groupName: string;
}

export default function XAxisChart({
  geneNames,
  leftOffset,
  groupName
}: Props): JSX.Element {
  const dispatch = useContext(DispatchContext);
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

  const tissue = groupName.split('--').at(0);
  const cellTypeName = groupName.split('--').at(-1);
  const headerName = `${cellTypeName} (${tissue})`;
  return (
    <XAxisWrapper
      width={heatmapWidth}
      left={leftOffset}    
    >
      {groupName && <GeneGroupWrapper width={heatmapWidth}>
        <MarkerGeneHeader>Marker genes</MarkerGeneHeader>
        <MarkerGeneHeaderButton
          onClick={() => dispatch?.(removeSelectedGenes(groupName))}
        >
          <Image
            src={xSvg.src}
            width="10"
            height="10"
            alt="remove marker gene group"
          />          
        </MarkerGeneHeaderButton>
        <GeneGroupName>{capitalize(headerName)}</GeneGroupName>
      </GeneGroupWrapper>  }  
      <XAxisContainer
        data-test-id="gene-labels"
        width={heatmapWidth}
        ref={xAxisRef}
      />
    </XAxisWrapper>
  );
}


function capitalize(str: string): string {
  return str.charAt(0).toUpperCase() + str.slice(1);
}
