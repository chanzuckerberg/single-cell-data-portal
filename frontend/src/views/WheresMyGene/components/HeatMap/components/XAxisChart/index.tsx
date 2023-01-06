import { useEffect, useState } from "react";
import { useDeleteGenesAndCellTypes } from "../../hooks/useDeleteGenesAndCellTypes";
import { CHART_LEFT_PADDING } from "../../style";
import { getHeatmapWidth, Y_AXIS_CHART_WIDTH_PX } from "../../utils";
import { XAxisWrapper, XAxisContainer, XAxisLabel } from "./style";

interface Props {
  geneNames: string[];
}

export default function XAxisChart({ geneNames }: Props): JSX.Element {
  const [heatmapWidth, setHeatmapWidth] = useState(getHeatmapWidth(geneNames));
  const { genesToDelete, handleGeneClick } = useDeleteGenesAndCellTypes();

  // Update heatmap size
  useEffect(() => {
    setHeatmapWidth(getHeatmapWidth(geneNames));
  }, [geneNames]);

  return (
    <XAxisWrapper
      width={heatmapWidth}
      left={Y_AXIS_CHART_WIDTH_PX + CHART_LEFT_PADDING}
    >
      <XAxisContainer data-test-id="gene-labels" width={heatmapWidth}>
        {geneNames.map((geneName) => (
          <XAxisLabel
            key={geneName}
            active={genesToDelete.includes(geneName)}
            data-test-id={`gene-label-${geneName}`}
            onClick={() => handleGeneClick(geneName)}
          >
            {geneName}
          </XAxisLabel>
        ))}
      </XAxisContainer>
    </XAxisWrapper>
  );
}
