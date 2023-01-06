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
  Y_AXIS_CHART_WIDTH_PX,
  formatLabel
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
import { SELECTED_STYLE } from "../../style";

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
  tissue,
}: {
  name: string;
  isMarkerGenes: boolean;
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
        {!FMG_EXCLUDE_TISSUES.includes(tissue) && (
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
          </InfoButtonWrapper>
        )}
      </FlexRow>
      <CellCountLabelStyle>{countString}</CellCountLabelStyle>
    </FlexRowJustified>
  );
};

function capitalize(str: string): string {
  return str.charAt(0).toUpperCase() + str.slice(1);
}