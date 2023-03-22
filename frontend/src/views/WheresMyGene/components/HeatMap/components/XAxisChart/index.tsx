import { Icon } from "czifui";
import React, { useContext, useEffect, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { DispatchContext } from "src/views/WheresMyGene/common/store";
import {
  deleteSingleGene,
  clearCellInfoCellType,
} from "src/views/WheresMyGene/common/store/actions";
import { useDeleteGenes } from "../../hooks/useDeleteGenes";
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
  InfoButtonWrapper,
} from "./style";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../../../GeneSearchBar/components/SaveExport";
import { StyledImage } from "../YAxisChart/style";
import InfoSVG from "../YAxisChart/icons/info-sign-icon.svg";
interface Props {
  geneNames: string[];
  generateGeneInfo: (gene: string) => void;
}

function GeneButton({
  geneName,
  genesToDelete,
  generateGeneInfo,
}: {
  geneName: string;
  genesToDelete: string[];
  handleGeneClick: (gene: string) => void;
  generateGeneInfo: (gene: string) => void;
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

        <InfoButtonWrapper
          className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}
          onClick={() => {
            generateGeneInfo(geneName);

            track(EVENTS.WMG_GENE_INFO, {
              gene: geneName,
            });

            if (!dispatch) return;

            // Clear cell type info here so that FMG panel closes
            dispatch(clearCellInfoCellType());
          }}
        >
          <StyledImage
            id="gene-info-button"
            src={InfoSVG.src}
            width="10"
            height="10"
            alt={`display gene info for ${geneName}`}
          />
        </InfoButtonWrapper>

      </XAxisLabel>
    </GeneButtonStyle>
  );
}

export default function XAxisChart({
  geneNames,
  generateGeneInfo,
}: Props): JSX.Element {
  const { genesToDelete, handleGeneClick } = useDeleteGenes();
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
            generateGeneInfo={generateGeneInfo}
          />
        ))}
      </XAxisContainer>
    </XAxisWrapper>
  );
}
