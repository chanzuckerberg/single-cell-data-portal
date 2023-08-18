import { Icon } from "@czi-sds/components";
import React, {
  Dispatch,
  SetStateAction,
  useContext,
  useEffect,
  useState,
} from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import {
  deleteSingleGene,
  selectGeneInfoFromXAxis,
  setXAxisHeight,
} from "src/views/WheresMyGene/common/store/actions";
import { useDeleteGenes } from "../../hooks/useDeleteGenes";
import { CHART_PADDING_PX, SELECTED_STYLE } from "../../style";
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
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../../../GeneSearchBar/components/SaveExport";
import { StyledImage } from "../YAxisChart/style";
import InfoSVG from "../YAxisChart/icons/info-sign-icon.svg";
import {
  X_AXIS_CHART_HEIGHT_PX,
  X_AXIS_HOVER_CONTAINER_HEIGHT_PX,
} from "src/views/WheresMyGene/common/constants";

interface Props {
  geneNames: string[];
}

export const GENE_LABEL_HOVER_CONTAINER_ID = "gene-hover-container";

function GeneButton({
  geneName,
  genesToDelete,
  currentMaxLabelHeight,
  setCurrentMaxLabelHeight,
}: {
  geneName: string;
  genesToDelete: string[];
  handleGeneClick: (gene: string) => void;
  currentMaxLabelHeight: number;
  setCurrentMaxLabelHeight: Dispatch<SetStateAction<number>>;
}): JSX.Element {
  const dispatch = useContext(DispatchContext);

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

  useEffect(() => {
    if (dispatch && length && length > currentMaxLabelHeight) {
      setCurrentMaxLabelHeight(length);
      dispatch(setXAxisHeight(length + X_AXIS_HOVER_CONTAINER_HEIGHT_PX)); // Account for hover icons
    }
  }, [dispatch, length, currentMaxLabelHeight, setCurrentMaxLabelHeight]);

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
          active={genesToDelete.includes(geneName)}
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

export default function XAxisChart({ geneNames }: Props): JSX.Element {
  const { genesToDelete, handleGeneClick } = useDeleteGenes();
  const [heatmapWidth, setHeatmapWidth] = useState(getHeatmapWidth(geneNames));
  const { xAxisHeight } = useContext(StateContext);

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
      <XAxisContainer
        data-testid="gene-labels"
        width={heatmapWidth}
        height={xAxisHeight}
      >
        {geneNames.map((geneName) => (
          <GeneButton
            key={geneName}
            geneName={geneName}
            genesToDelete={genesToDelete}
            handleGeneClick={handleGeneClick}
            currentMaxLabelHeight={currentMaxLabelHeight}
            setCurrentMaxLabelHeight={setCurrentMaxLabelHeight}
          />
        ))}
      </XAxisContainer>
    </XAxisWrapper>
  );
}
