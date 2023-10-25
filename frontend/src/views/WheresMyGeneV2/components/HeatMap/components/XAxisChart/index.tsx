import { Icon } from "@czi-sds/components";
import React, { useContext, useEffect, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGeneV2/common/store";
import {
  deleteSingleGene,
  selectGeneInfoFromXAxis,
  setXAxisHeight,
} from "src/views/WheresMyGeneV2/common/store/actions";
import { useDeleteGenes } from "../../hooks/useDeleteGenes";
import { CHART_PADDING_PX, SELECTED_STYLE } from "../../style_old";
import {
  getHeatmapWidth,
  formatLabel,
  Y_AXIS_CHART_WIDTH_PX,
} from "../../utils";
import {
  XAxisWrapper,
  XAxisContainer,
  XAxisLabel,
  GeneButtonStyle,
  XAxisGeneName,
  InfoButtonWrapper,
  HoverContainer,
  DeleteButtonWrapper,
} from "./style";
import { StyledImage } from "../YAxisChart/style";
import InfoSVG from "src/common/images/info-sign-icon.svg";
import {
  X_AXIS_CHART_HEIGHT_PX,
  X_AXIS_HOVER_CONTAINER_HEIGHT_PX,
} from "src/views/WheresMyGeneV2/common/constants";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "src/views/WheresMyGeneV2/components/GeneSearchBar/components/SaveExport";
import GeneSearchBar from "src/views/WheresMyGeneV2/components/GeneSearchBar";

interface Props {
  geneNames: string[];
  sidebarWidth: number;
}

export const GENE_LABEL_HOVER_CONTAINER_ID = "gene-hover-container";

function GeneButton({
  geneName,
  active,
  currentFont,
  formattedLabel,
}: {
  active: boolean;
  formattedLabel: string;
  currentFont: string;
  geneName: string;
  handleGeneClick: (gene: string) => void;
}): JSX.Element {
  const dispatch = useContext(DispatchContext);

  return (
    <GeneButtonStyle
      id={`gene-label-${geneName}`}
      data-testid={`gene-label-${geneName}`}
    >
      <HoverContainer
        id={GENE_LABEL_HOVER_CONTAINER_ID}
        data-testid={GENE_LABEL_HOVER_CONTAINER_ID}
      >
        <DeleteButtonWrapper
          data-testid={`gene-delete-icon-${geneName}`}
          className={`gene-delete-icon ${EXCLUDE_IN_SCREENSHOT_CLASS_NAME}`}
          onClick={() => {
            track(EVENTS.WMG_DELETE_GENE, { gene: geneName });

            if (dispatch) {
              dispatch(deleteSingleGene(geneName));
            }
          }}
        >
          <Icon sdsIcon="trashCan" sdsSize="s" sdsType="interactive"></Icon>
        </DeleteButtonWrapper>

        <InfoButtonWrapper
          data-testid="gene-info-button-x-axis"
          className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}
          onClick={() => {
            if (!dispatch) return;
            track(EVENTS.WMG_GENE_INFO, {
              gene: geneName,
            });
            // This will populate gene info and clear cell type info so that FMG panel closes
            dispatch(selectGeneInfoFromXAxis(geneName));
          }}
        >
          <StyledImage
            src={InfoSVG.src}
            width="10"
            height="10"
            alt={`display gene info for ${geneName}`}
            data-testid={`gene-info-icon-${geneName}`}
          />
        </InfoButtonWrapper>
      </HoverContainer>
      <XAxisLabel className={"gene-label-container"}>
        <XAxisGeneName
          active={active}
          font={currentFont}
          className={`gene-name-${geneName}`}
          data-testid={`gene-name-${geneName}`}
        >
          {formattedLabel}
        </XAxisGeneName>
      </XAxisLabel>
    </GeneButtonStyle>
  );
}

export default function XAxisChart({
  geneNames,
  sidebarWidth,
}: Props): JSX.Element {
  const { genesToDelete, handleGeneClick } = useDeleteGenes();
  const [heatmapWidth, setHeatmapWidth] = useState(getHeatmapWidth(geneNames));
  const { xAxisHeight } = useContext(StateContext);
  const dispatch = useContext(DispatchContext);

  // This is used to calculate the current longest gene name
  const [currentMaxLabelHeight, setCurrentMaxLabelHeight] = useState(0);

  // Update heatmap size
  useEffect(() => {
    setHeatmapWidth(getHeatmapWidth(geneNames));
    setCurrentMaxLabelHeight(0); // Reset when selected gene labels change
  }, [geneNames]);

  return (
    <XAxisWrapper
      width={heatmapWidth}
      left={Y_AXIS_CHART_WIDTH_PX + CHART_PADDING_PX}
      height={xAxisHeight}
    >
      <GeneSearchBar
        sidebarWidth={sidebarWidth}
        className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}
      />
      <XAxisContainer
        data-testid="gene-labels"
        width={heatmapWidth}
        height={xAxisHeight}
      >
        {geneNames.map((geneName) => {
          const active = genesToDelete.includes(geneName);
          const currentFont = `
          normal
          ${active ? SELECTED_STYLE.fontWeight : "normal"}
          ${SELECTED_STYLE.fontSize}px ${SELECTED_STYLE.fontFamily}
        `;
          const { text: formattedLabel, length } = formatLabel(
            geneName,
            X_AXIS_CHART_HEIGHT_PX, // Gene label length is capped to this value
            currentFont
          );

          if (length > currentMaxLabelHeight) {
            setCurrentMaxLabelHeight(length);
            if (dispatch) {
              dispatch(
                setXAxisHeight(length + X_AXIS_HOVER_CONTAINER_HEIGHT_PX)
              );
            }
          }

          return (
            <GeneButton
              key={geneName}
              currentFont={currentFont}
              formattedLabel={formattedLabel}
              geneName={geneName}
              active={active}
              handleGeneClick={handleGeneClick}
            />
          );
        })}
      </XAxisContainer>
    </XAxisWrapper>
  );
}
