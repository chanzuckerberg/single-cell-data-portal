import { Icon } from "czifui";
import React, { useContext, useEffect, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { DispatchContext } from "src/views/WheresMyGene/common/store";
import { deleteSingleGene } from "src/views/WheresMyGene/common/store/actions";
import { useDeleteGenesAndCellTypes } from "../../hooks/useDeleteGenesAndCellTypes";
import { CHART_LEFT_PADDING, SELECTED_STYLE } from "../../style";
import {
  getHeatmapWidth,
  formatLabel,
  X_AXIS_CHART_HEIGHT_PX,
  Y_AXIS_CHART_WIDTH_PX,
} from "../../utils";
import {
  XAxisWrapper,
  XAxisContainer,
  XAxisLabel,
  GeneButtonStyle,
  XAxisGeneName,
} from "./style";

interface Props {
  geneNames: string[];
}

function GeneButton({
  geneName,
  genesToDelete,
}: {
  geneName: string;
  genesToDelete: string[];
  handleGeneClick: (gene: string) => void;
}): JSX.Element {
  const dispatch = useContext(DispatchContext);

  const active = genesToDelete.includes(geneName);
  const currentFont = `
    normal
    ${active ? SELECTED_STYLE.fontWeight : "normal"}
    ${SELECTED_STYLE.fontSize}px ${SELECTED_STYLE.fontFamily}
  `;
  const formattedLabel = formatLabel(
    geneName,
    X_AXIS_CHART_HEIGHT_PX,
    currentFont
  );

  return (
    <GeneButtonStyle data-test-id={`gene-label-${geneName}`}>
      <XAxisLabel className={"gene-label-container"}>
        <div
          data-test-id={"gene-delete-button"}
          className="gene-delete-icon"
          onClick={() => {
            track(EVENTS.WMG_DELETE_GENE, { gene: geneName });

            if (dispatch) {
              dispatch(deleteSingleGene(geneName));
            }
          }}
        >
          <Icon sdsIcon="trashCan" sdsSize="s" sdsType="interactive"></Icon>
        </div>
        <XAxisGeneName
          active={genesToDelete.includes(geneName)}
          font={currentFont}
        >
          {formattedLabel}
        </XAxisGeneName>
      </XAxisLabel>
    </GeneButtonStyle>
  );
}

export default function XAxisChart({ geneNames }: Props): JSX.Element {
  const { genesToDelete, handleGeneClick } = useDeleteGenesAndCellTypes();
  const [heatmapWidth, setHeatmapWidth] = useState(getHeatmapWidth(geneNames));

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
          <GeneButton
            key={geneName}
            geneName={geneName}
            genesToDelete={genesToDelete}
            handleGeneClick={handleGeneClick}
          />
        ))}
      </XAxisContainer>
    </XAxisWrapper>
  );
}
