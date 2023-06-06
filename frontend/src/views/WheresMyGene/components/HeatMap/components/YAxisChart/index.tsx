import { memo, useEffect, useMemo, useState } from "react";
import { CellType, Tissue } from "src/views/WheresMyGene/common/types";
import {
  CellTypeMetadata,
  deserializeCellTypeMetadata,
  formatLabel,
  getAllSerializedCellTypeMetadata,
  getHeatmapHeight,
  hyphenize,
  Y_AXIS_CHART_WIDTH_PX,
} from "../../utils";
import InfoSVG from "./icons/info-sign-icon.svg";
import {
  CellCountLabelStyle,
  CellTypeLabelStyle,
  CellTypeLabelTooltipStyle,
  Container,
  FlexRow,
  FlexRowJustified,
  HiddenCellTypeLabelStyle,
  StyledImage,
  TissueName,
  TissueWrapper,
  Wrapper,
} from "./style";
import { SELECTED_STYLE } from "../../style";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../../../GeneSearchBar/components/SaveExport";
import { COMPARE_OPTION_ID_FOR_AGGREGATED } from "src/common/queries/wheresMyGene";
import { InfoButtonWrapper } from "src/components/common/Filter/common/style";
import { Tooltip } from "@czi-sds/components";

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
  const tissueKey = hyphenize(tissue);

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
    <Wrapper>
      <TissueWrapper height={heatmapHeight}>
        <TissueName>{capitalize(tissue)}</TissueName>
      </TissueWrapper>
      <Container
        data-testid={`cell-type-labels-${tissueKey}`}
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

            const { text: paddedName } = formatLabel(
              name,
              Y_AXIS_CHART_WIDTH_PX - 90, // scale based on y-axis width
              selectedFont // prevents selected style from overlapping count
            );
            return (
              <CellTypeButton
                key={`${cellType}-cell-type-button`}
                formattedName={paddedName}
                name={name}
                metadata={cellType}
                tissueID={tissueID}
                tissue={tissue}
                generateMarkerGenes={generateMarkerGenes}
                data-testid="cell-type-label"
              />
            );
          })}
      </Container>
    </Wrapper>
  );
});

const CellTypeButton = ({
  formattedName,
  name,
  metadata,
  generateMarkerGenes,
  tissueID,
  tissue,
}: {
  formattedName: string;
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

  const isTruncated = formattedName.includes("...");

  return (
    <FlexRowJustified
      id="cell-type-label-count"
      data-testid="cell-type-label-count"
    >
      <FlexRow>
        <CellTypeLabelStyle>
          <Tooltip
            title={
              // Set tooltip content only if name is truncated
              isTruncated ? (
                <CellTypeLabelTooltipStyle data-testid="cell-type-name-tooltip">
                  {name}
                </CellTypeLabelTooltipStyle>
              ) : null
            }
            sdsStyle="light"
            arrow
            placement="left"
            enterNextDelay={700}
          >
            {/* Must be wrapped in div and not fragment or else tooltip content won't render */}
            <div>
              {/* Hidden labels are only needed if name is truncated */}
              {isTruncated && (
                <HiddenCellTypeLabelStyle data-testid="cell-type-full-name">
                  {name}
                </HiddenCellTypeLabelStyle>
              )}
              <div id="cell-type-name" data-testid="cell-type-name">
                {formattedName}
              </div>
            </div>
          </Tooltip>
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
                data-testid="marker-gene-button"
                src={InfoSVG.src}
                width="10"
                height="10"
                alt={`display marker genes for ${cellType.name}`}
              />
            </InfoButtonWrapper>
          )}
      </FlexRow>
      <CellCountLabelStyle id="cell-count" data-testid="cell-count">
        {countString}
      </CellCountLabelStyle>
    </FlexRowJustified>
  );
};

export function capitalize(str: string): string {
  return str.charAt(0).toUpperCase() + str.slice(1);
}
