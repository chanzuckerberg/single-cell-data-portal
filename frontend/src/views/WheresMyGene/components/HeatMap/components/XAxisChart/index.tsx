import { Icon } from "czifui";
import React, { useContext, useEffect, useState } from "react";
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
} from "./style";

interface Props {
  geneNames: string[];
}

function GeneButton({
  geneName,
  genesToDelete,
  handleGeneClick,
}: {
  geneName: string;
  genesToDelete: string[];
  handleGeneClick: (gene: string) => void;
}): JSX.Element {
  const dispatch = useContext(DispatchContext);

  const [hideTrash, setHideTrashClass] = useState("hide-gene-delete");

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
    <GeneButtonStyle
      onClick={() => handleGeneClick(geneName)}
      data-test-id={`gene-label-${geneName}`}
      onMouseEnter={() => {
        setHideTrashClass("");
      }}
      onMouseLeave={() => {
        setHideTrashClass("hide-gene-delete");
      }}
    >
      <XAxisLabel active={genesToDelete.includes(geneName)} font={currentFont}>
        <div
          data-test-id={"gene-delete-button"}
          onClick={() => {
            if (dispatch) {
              dispatch(deleteSingleGene(geneName));
            }
          }}
        >
          <Icon
            sdsIcon="trashCan"
            sdsSize="s"
            sdsType="interactive"
            className={hideTrash}
          ></Icon>
        </div>
        <div>{formattedLabel}</div>
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
