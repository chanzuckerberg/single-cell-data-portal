import { memo, useContext, useEffect, useMemo, useState } from "react";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import Image from "next/image";
import { DispatchContext } from "src/views/WheresMyGene/common/store";
import { resetTissueCellTypes } from "src/views/WheresMyGene/common/store/actions";
import { CellType, Tissue } from "src/views/WheresMyGene/common/types";
import { useDeleteGenesAndCellTypes } from "../../hooks/useDeleteGenesAndCellTypes";
import {
  CellTypeMetadata,
  deserializeCellTypeMetadata,
  getAllSerializedCellTypeMetadata,
  getHeatmapHeight,
  SELECTED_STYLE,
  Y_AXIS_CHART_WIDTH_PX,
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
  CellTypeButtonStyle,
  CellCountLabelStyle,
  FlexRowJustified,
  FlexRow,
  InfoButtonWrapper,
} from "./style";

const MAX_DEPTH = 2;

interface Props {
  cellTypes?: CellType[];
  hasDeletedCellTypes: boolean;
  availableCellTypes: CellType[];
  tissue: Tissue;
  tissueID: string;
  generateMarkerGenes: (cellType: CellType, tissueID: string) => void;
  selectedOrganismId: string;
}

// List of Tissues to exclude from FMG
const FMG_EXCLUDE_TISSUES = [
  "blood"
];

type Coord = [number, number];
export default memo(function YAxisChart({
  cellTypes = [],
  hasDeletedCellTypes,
  availableCellTypes,
  tissue,
  generateMarkerGenes,
  tissueID,
}: Props): JSX.Element {
  const tissueKey = tissue.replace(/\s+/g, "-");

  const dispatch = useContext(DispatchContext);
  const isMarkerGenes = get(FEATURES.MARKER_GENES) === BOOLEAN.TRUE;

  const { handleCellTypeClick } = useDeleteGenesAndCellTypes();

  const [heatmapHeight, setHeatmapHeight] = useState(
    getHeatmapHeight(cellTypes)
  );

  // Update heatmap size
  useEffect(() => {
    setHeatmapHeight(getHeatmapHeight(cellTypes));
  }, [cellTypes]);

  const cellTypeMetadata = useMemo(() => {
    return getAllSerializedCellTypeMetadata(cellTypes, tissue);
  }, [cellTypes, tissue]);

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
      >
        {cellTypeMetadata
          .slice()
          .reverse()
          .map((cellType) => {
            const { name, depth = 0 } = deserializeCellTypeMetadata(
              cellType as CellTypeMetadata
            );
            const displayDepth = Math.min(depth, MAX_DEPTH);

            const { fontWeight, fontSize, fontFamily } = SELECTED_STYLE;
            const selectedFont = `${fontWeight} ${fontSize}px ${fontFamily}`;

            const paddedName = formatCellLabel(
              name,
              Y_AXIS_CHART_WIDTH_PX - 90, // scale based on y-axis width
              selectedFont, // prevents selected style from overlapping count
              displayDepth
            );
            return (
              <CellTypeButton
                key={`${cellType}-cell-type-button`}
                name={paddedName}
                metadata={cellType}
                onClick={() => handleCellTypeClick(cellType)}
                tissueID={tissueID}
                generateMarkerGenes={generateMarkerGenes}
                isMarkerGenes={isMarkerGenes}
              />
            );
          })}
      </Container>
    </Wrapper>
  );

  function handleResetTissue(tissue: Tissue) {
    if (!dispatch) return;

    dispatch(resetTissueCellTypes(tissue, availableCellTypes));
  }
});

const CellTypeButton = ({
  name,
  metadata,
  onClick,
  isMarkerGenes,
  generateMarkerGenes,
  tissueID,
}: {
  name: string;
  isMarkerGenes: boolean;
  metadata: CellTypeMetadata;
  onClick: () => void;
  generateMarkerGenes: (cellType: CellType, tissueID: string) => void;
  tissueID: string;
}) => {
  const [active, setActive] = useState(false);
  useEffect(() => {
    setActive(false);
  }, [metadata]);

  const { total_count } = deserializeCellTypeMetadata(metadata);
  const formattedString = Intl.NumberFormat("en-US", {
    maximumFractionDigits: 1,
    notation: "compact",
  }).format(total_count);
  const countString = `${formattedString}${
    formattedString !== total_count.toString() ? "+" : ""
  }`;

  return (
    <FlexRowJustified>
      <FlexRow>
        <CellTypeButtonStyle
          active={active}
          onClick={() => {
            setActive(!active);
            onClick();
          }}
        >
          {name}
        </CellTypeButtonStyle>
        {!FMG_EXCLUDE_TISSUES.includes(tissue) &&
         <InfoButtonWrapper
          style={{
            paddingTop: "3px",
            cursor: "pointer",
          }}
          onClick={() => {
            if (isMarkerGenes) {
              const cellType = deserializeCellTypeMetadata(metadata);
              generateMarkerGenes(cellType, tissueID);
            }
          }}
        >
          <StyledImage
            id={"marker-gene-button"}
            src={InfoSVG.src}
            width="10"
            height="10"
            alt={`display marker genes for ${
              deserializeCellTypeMetadata(metadata).name
            }`}
          />
        </InfoButtonWrapper>}
      </FlexRow>
      <CellCountLabelStyle>{countString}</CellCountLabelStyle>
    </FlexRowJustified>
  );
};

function capitalize(str: string): string {
  return str.charAt(0).toUpperCase() + str.slice(1);
}

/**
 * Used to calculate text pixel widths. Should be only created once.
 */
const CTX =
  (typeof document !== "undefined" &&
    document.createElement("canvas").getContext("2d")) ||
  null;

/**
 * Formats and truncates the cell type name to a given width
 *
 * @param name The text to truncate
 * @param maxWidth The max width in pixels the string should be
 * @param font The font family and font size as a string. Ex. "bold 12px sans-serif"
 * @param displayDepth The depth of the cell type name (indentation/padding)
 * @returns The string fixed to a certain pixel width
 */
function formatCellLabel(
  name: string,
  maxWidth: number,
  font: string,
  displayDepth = 0
): string {
  CTX!.font = font;
  const ellipsisWidth = CTX!.measureText("...").width;

  const padding = " ".repeat(displayDepth * 8);
  const paddingWidth = CTX!.measureText(padding).width;

  if (CTX!.measureText(name).width + paddingWidth <= maxWidth) {
    return padding + name;
  }

  const labelHalfWidth = (maxWidth - paddingWidth - ellipsisWidth) / 2;

  const firstHalf = getFixedWidth(name, labelHalfWidth, font);
  const secondHalf = getFixedWidth(name, labelHalfWidth, font, true);

  return padding + firstHalf + "..." + secondHalf;
}

/**
 * Truncates the string to a given width
 *
 * @param text The text to truncate
 * @param maxWidth The max width in pixels the string should be
 * @param font The font family and font size as a string. Ex. "bold 12px sans-serif"
 * @param reverse Whether to truncate the end or beginning of the string
 * @returns The string fixed to a certain pixel width
 */
export function getFixedWidth(
  text: string,
  maxWidth: number,
  font: string,
  reverse = false
): string {
  CTX!.font = font;

  if (reverse) {
    for (let i = text.length; i >= 0; i--) {
      const substring = text.substring(i - 1);
      const textWidth = CTX!.measureText(substring).width;
      if (textWidth > maxWidth) {
        return text.substring(i);
      }
    }
  } else {
    for (let i = 0; i < text.length; i++) {
      const substring = text.substring(0, i + 1);
      const textWidth = CTX!.measureText(substring).width;
      if (textWidth > maxWidth) {
        return text.substring(0, i);
      }
    }
  }

  return text;
}
