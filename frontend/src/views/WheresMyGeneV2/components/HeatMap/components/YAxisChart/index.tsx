import { memo, useCallback, useEffect, useMemo, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../../../GeneSearchBar/components/SaveExport";
import { InfoButtonWrapper } from "src/components/common/Filter/common/style";
import { Icon, Tooltip } from "@czi-sds/components";
import InfoSVG from "src/views/WheresMyGene/components/HeatMap/components/YAxisChart/icons/info-sign-icon.svg";

import {
  COMPARE_OPTION_ID_FOR_AGGREGATED,
  CellTypeRow,
} from "src/common/queries/wheresMyGene";
import { CellType, Tissue } from "src/views/WheresMyGene/common/types";
import {
  CellTypeMetadata,
  Y_AXIS_CHART_WIDTH_PX,
  deserializeCellTypeMetadata,
  formatLabel,
  getAllSerializedCellTypeMetadata,
  getHeatmapHeight,
  hyphenize,
} from "src/views/WheresMyGene/components/HeatMap/utils";
import { SELECTED_STYLE } from "src/views/WheresMyGene/components/HeatMap/style";
import {
  CellCountLabelStyle,
  CellTypeLabelStyle,
  CellTypeLabelTooltipStyle,
  Container,
  FlexRow,
  FlexRowJustified,
  HiddenCellTypeLabelStyle,
  StyledImage,
  TissueHeaderLabelStyle,
  Wrapper,
  TissueLabel,
} from "src/views/WheresMyGene/components/HeatMap/components/YAxisChart/style";

interface Props {
  cellTypes: CellTypeRow[];
  tissue: Tissue;
  tissueID: string;
  generateMarkerGenes: (cellType: CellType, tissueID: string) => void;
  handleExpandCollapse: (tissueID: string, tissueName: Tissue) => void;
  expandedTissues: Set<Tissue>;
  selectedOrganismId: string;
}

// List of Tissues to exclude from FMG
const FMG_EXCLUDE_TISSUES = ["blood"];

export default memo(function YAxisChart({
  cellTypes = [],
  tissue,
  generateMarkerGenes,
  handleExpandCollapse,
  expandedTissues,
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
    <Wrapper id={`${hyphenize(tissue)}-y-axis`}>
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
            const expanded = expandedTissues.has(tissueID);

            const { text: paddedName } = formatLabel(
              name,
              Y_AXIS_CHART_WIDTH_PX - 90, // scale based on y-axis width
              selectedFont // prevents selected style from overlapping count
            );
            if (name === tissue) {
              return (
                <TissueHeaderButton
                  key={`${cellType}-cell-type-button`}
                  formattedName={paddedName}
                  metadata={cellType}
                  tissueID={tissueID}
                  generateMarkerGenes={generateMarkerGenes}
                  expanded={expanded}
                  tissueName={name}
                  handleExpandCollapse={handleExpandCollapse}
                  data-testid="cell-type-label"
                />
              );
            }
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

const TissueHeaderButton = ({
  formattedName,
  metadata,
  handleExpandCollapse,
  expanded,
  tissueID,
  tissueName,
}: {
  formattedName: string;
  metadata: CellTypeMetadata;
  generateMarkerGenes: (cellType: CellType, tissueID: string) => void;
  tissueID: string;
  tissueName: Tissue;
  handleExpandCollapse: (tissueID: string, tissueName: Tissue) => void;
  expanded: boolean;
}) => {
  const { total_count } = deserializeCellTypeMetadata(metadata);
  const formattedString = Intl.NumberFormat("en-US", {
    maximumFractionDigits: 1,
    notation: "compact",
  }).format(total_count);
  const countString = `${formattedString}`.toLowerCase();

  return (
    <FlexRowJustified
      className="cell-type-label-count"
      data-testid="cell-type-label-count"
    >
      <FlexRow
        onClick={useCallback(() => {
          handleExpandCollapse(tissueID, tissueName);
        }, [tissueID, tissueName, handleExpandCollapse])}
      >
        <Icon
          sdsIcon={expanded ? "triangleDown" : "triangleRight"}
          sdsSize="s"
          color="gray"
          sdsType="static"
        />
        <TissueHeaderLabelStyle>
          <div>
            <TissueLabel
              className="cell-type-name"
              data-testid="cell-type-name"
            >
              {capitalize(formattedName)}
            </TissueLabel>
          </div>
        </TissueHeaderLabelStyle>
      </FlexRow>
      <CellCountLabelStyle className="cell-count" data-testid="cell-count">
        {countString}
      </CellCountLabelStyle>
    </FlexRowJustified>
  );
};
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
      className="cell-type-label-count"
      data-testid="cell-type-label-count"
    >
      <FlexRow>
        <CellTypeLabelStyle>
          <Tooltip
            title={
              // Set tooltip content only if name is truncated
              isTruncated ? (
                <CellTypeLabelTooltipStyle>{name}</CellTypeLabelTooltipStyle>
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
              <div className="cell-type-name" data-testid="cell-type-name">
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
                /**
                 * (thuang): https://nextjs.org/docs/pages/api-reference/components/image-legacy#layout
                 * Use the <StyledImage /> width and height, since default is `intrinsic`
                 */
                layout="fixed"
                width="10"
                height="10"
                alt={`display marker genes for ${cellType.name}`}
              />
            </InfoButtonWrapper>
          )}
      </FlexRow>
      <CellCountLabelStyle className="cell-count" data-testid="cell-count">
        {countString}
      </CellCountLabelStyle>
    </FlexRowJustified>
  );
};

export function capitalize(str: string): string {
  return str.charAt(0).toUpperCase() + str.slice(1);
}
