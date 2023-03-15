import { memo, useEffect, useMemo, useState } from "react";
import { CellType, Tissue } from "src/views/WheresMyGene/common/types";
import {
  CellTypeMetadata,
  deserializeCellTypeMetadata,
  formatLabel,
  getAllSerializedCellTypeMetadata,
  getHeatmapHeight,
  Y_AXIS_CHART_WIDTH_PX,
} from "../../utils";
import InfoSVG from "./icons/info-sign-icon.svg";
import {
  CellCountLabelStyle,
  CellTypeLabelStyle,
  Container,
  FlexRow,
  FlexRowJustified,
  InfoButtonWrapper,
  StyledImage,
  TissueName,
  TissueWrapper,
  Wrapper,
} from "./style";
import { SELECTED_STYLE } from "../../style";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../../../GeneSearchBar/components/SaveImage";
import { COMPARE_OPTION_ID_FOR_AGGREGATED } from "src/common/queries/wheresMyGene";

interface Props {
  cellTypes?: CellType[];
  tissue: Tissue;
  tissueID: string;
  generateMarkerGenes: (cellType: CellType, tissueID: string) => void;
  selectedOrganismId: string;
}

// List of Tissues to exclude from FMG
const FMG_EXCLUDE_TISSUES = ["blood"];

export default memo(function YAxisChart({
  cellTypes = [],
  tissue,
  generateMarkerGenes,
  tissueID,
}: Props): JSX.Element {
  const tissueKey = tissue.replace(/\s+/g, "-");

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
    <Wrapper id={`${tissue.replace(/\s+/g, "-")}-y-axis`}>
      <TissueWrapper height={heatmapHeight}>
        <TissueName>{capitalize(tissue)}</TissueName>
      </TissueWrapper>
      <Container
        data-test-id={`cell-type-labels-${tissueKey}`}
        height={heatmapHeight}
      >
        {cellTypeMetadata
          .slice()
          .reverse()
          .map((cellType) => {
            const { name } = deserializeCellTypeMetadata(
              cellType as CellTypeMetadata
            );
            const { fontWeight, fontSize, fontFamily } = SELECTED_STYLE;
            const selectedFont = `${fontWeight} ${fontSize}px ${fontFamily}`;

            const paddedName = formatLabel(
              name,
              Y_AXIS_CHART_WIDTH_PX - 90, // scale based on y-axis width
              selectedFont // prevents selected style from overlapping count
            );
            return (
              <CellTypeButton
                key={`${cellType}-cell-type-button`}
                name={paddedName}
                metadata={cellType}
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
});

const CellTypeButton = ({
  name,
  metadata,
  generateMarkerGenes,
  tissueID,
  tissue,
}: {
  name: string;
  metadata: CellTypeMetadata;
  generateMarkerGenes: (cellType: CellType, tissueID: string) => void;
  tissueID: string;
  tissue: Tissue;
}) => {
  const { total_count } = deserializeCellTypeMetadata(metadata);
  const formattedString = Intl.NumberFormat("en-US", {
    maximumFractionDigits: 1,
    notation: "compact",
  }).format(total_count);
  const countString = `${formattedString}`.toLowerCase();

  const cellType = deserializeCellTypeMetadata(metadata);

  return (
    <FlexRowJustified data-test-id="cell-type-label-count">
      <FlexRow>
        <CellTypeLabelStyle data-test-id="cell-type-name">
          {name}
        </CellTypeLabelStyle>

        {!FMG_EXCLUDE_TISSUES.includes(tissue) &&
          cellType &&
          cellType.total_count > 25 &&
          cellType.optionId === COMPARE_OPTION_ID_FOR_AGGREGATED && (
            <InfoButtonWrapper
              className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}
              style={{
                cursor: "pointer",
                margin: "auto",
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
                id="marker-gene-button"
                src={InfoSVG.src}
                width="10"
                height="10"
                alt={`display marker genes for ${cellType.name}`}
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
