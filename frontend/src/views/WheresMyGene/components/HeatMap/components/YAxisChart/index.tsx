import Image from "next/image";
import { memo, useContext, useEffect, useMemo, useState } from "react";
import { DispatchContext } from "src/views/WheresMyGene/common/store";
import { resetTissueCellTypes } from "src/views/WheresMyGene/common/store/actions";
import { CellType, Tissue } from "src/views/WheresMyGene/common/types";
import { useDeleteGenesAndCellTypes } from "../../hooks/useDeleteGenesAndCellTypes";
import {
  CellTypeMetadata,
  deserializeCellTypeMetadata,
  formatLabel,
  getAllSerializedCellTypeMetadata,
  getHeatmapHeight,
  Y_AXIS_CHART_WIDTH_PX,
} from "../../utils";
import InfoSVG from "./icons/info-sign-icon.svg";
import ReplaySVG from "./icons/replay.svg";
import {
  CellCountLabelStyle,
  CellTypeButtonStyle,
  Container,
  FlexRow,
  FlexRowJustified,
  InfoButtonWrapper,
  ResetImageWrapper,
  StyledImage,
  TissueName,
  TissueWrapper,
  Wrapper,
} from "./style";
import { SELECTED_STYLE } from "../../style";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../../../GeneSearchBar/components/SaveImage";

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
const FMG_EXCLUDE_TISSUES = ["blood"];

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
  
  const isRollup = get(FEATURES.IS_ROLLUP) === BOOLEAN.TRUE;
  return (
    <Wrapper id={`${tissue.replace(/\s+/g, "-")}-y-axis`}>
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
            const displayDepth = isRollup ? 0 : Math.min(depth, MAX_DEPTH);

            const { fontWeight, fontSize, fontFamily } = SELECTED_STYLE;
            const selectedFont = `${fontWeight} ${fontSize}px ${fontFamily}`;

            const paddedName = formatLabel(
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
                tissue={tissue}
                generateMarkerGenes={generateMarkerGenes}
                date-test-id="cell-type-label"
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
  generateMarkerGenes,
  tissueID,
  tissue,
}: {
  name: string;
  metadata: CellTypeMetadata;
  onClick: () => void;
  generateMarkerGenes: (cellType: CellType, tissueID: string) => void;
  tissueID: string;
  tissue: Tissue;
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

  const cellType = deserializeCellTypeMetadata(metadata);

  return (
    <FlexRowJustified data-test-id="cell-type-label-count">
      <FlexRow>
        <CellTypeButtonStyle
          active={active}
          onClick={() => {
            setActive(!active);
            onClick();
          }}
          data-test-id="cell-type-label"
        >
          {name}
        </CellTypeButtonStyle>
        {!FMG_EXCLUDE_TISSUES.includes(tissue) &&
          cellType &&
          cellType.total_count > 25 && (
            <InfoButtonWrapper
              className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}
              style={{
                cursor: "pointer",
                paddingTop: "3px",
              }}
              onClick={() => {
                if (cellType) {
                  generateMarkerGenes(cellType, tissueID);
                  track(EVENTS.WMG_FMG_INFO_CLICKED, {
                    combination: `${cellType.name}, ${tissue}}`,
                  });
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
            </InfoButtonWrapper>
          )}
      </FlexRow>
      <CellCountLabelStyle data-test-id="cell-count">
        {countString}
      </CellCountLabelStyle>
    </FlexRowJustified>
  );
};

function capitalize(str: string): string {
  return str.charAt(0).toUpperCase() + str.slice(1);
}
