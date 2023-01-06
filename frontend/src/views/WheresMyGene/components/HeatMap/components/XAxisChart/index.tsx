import { useEffect, useState } from "react";
import { useDeleteGenesAndCellTypes } from "../../hooks/useDeleteGenesAndCellTypes";
import { CHART_LEFT_PADDING, SELECTED_STYLE } from "../../style";
import {
  getHeatmapWidth,
  X_AXIS_CHART_HEIGHT_PX,
  Y_AXIS_CHART_WIDTH_PX,
} from "../../utils";
import {
  XAxisWrapper,
  XAxisContainer,
  XAxisLabel,
  GeneButtonStyle,
} from "./style";
import { formatLabel } from "../YAxisChart";

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
        {geneNames.map((geneName) => {
          const active = genesToDelete.includes(geneName);
          const currentFont = `
            normal
            ${active ? SELECTED_STYLE.fontWeight : "normal"}
            ${SELECTED_STYLE.fontSize}px ${SELECTED_STYLE.fontFamily}
          `;
          const formattedLabel = formatLabel(
            geneName,
            X_AXIS_CHART_HEIGHT_PX,
            currentFont,
            0
          );
          return (
            <GeneButtonStyle
              key={geneName}
              onClick={() => handleGeneClick(geneName)}
              data-test-id={`gene-label-${geneName}`}
            >
              <XAxisLabel
                active={genesToDelete.includes(geneName)}
                font={currentFont}
              >
                {formattedLabel}
              </XAxisLabel>
            </GeneButtonStyle>
          );
        })}
      </XAxisContainer>
    </XAxisWrapper>
  );
}
