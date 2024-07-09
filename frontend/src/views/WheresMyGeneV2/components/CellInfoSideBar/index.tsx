import { Icon, Tooltip } from "@czi-sds/components";
import React from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { ROUTES } from "src/common/constants/routes";
import { BetaChip } from "src/components/Header/style";
import {
  ButtonContainer,
  CopyGenesButton,
  NoMarkerGenesContainer,
  NoMarkerGenesDescription,
  NoMarkerGenesHeader,
  StyledMarkerGeneHeader,
  StyledTooltip,
  TooltipContent,
  TissueName,
  TooltipButton,
  ButtonWrapper,
  TooltipLink,
  AddToDotplotButton,
} from "./style";
import { StyledLink } from "src/components/GeneInfoSideBar/style";
import InfoSVG from "src/common/images/info-sign-icon.svg";
import { InfoButtonWrapper } from "src/components/common/Filter/common/style";
import { CellInfoBarProps } from "./types";
import {
  MARKER_GENES_TOOLTIP_CONTENT,
  MARKER_GENE_LABEL,
  MARKER_SCORE_HIGH_CONTENT,
  MARKER_SCORE_LABEL,
  MARKER_SCORE_LOW_CONTENT,
  MARKER_SCORE_MEDIUM_CONTENT,
  MARKER_SCORE_TOOLTIP_CONTENT,
  MARKER_SCORE_TOOLTIP_LINK_TEXT,
  EFFECT_SIZE,
  MARKER_SCORE_TOOLTIP_TEST_ID,
  SPECIFICITY_TOOLTIP_TEST_ID,
} from "src/common/constants/markerGenes";
import {
  TISSUES_WITHOUT_MARKER_GENES,
  MARKER_SCORE_CELLGUIDE_LINK_TEXT,
  MARKER_SCORE_DOTPLOT_BUTTON_TEXT,
  NO_MARKER_GENES_DESCRIPTION,
  NO_MARKER_GENES_FOR_BLOOD_DESCRIPTION,
  NO_MARKER_GENES_HEADER,
  TABLE_HEADER_GENE,
  TOO_FEW_CELLS_NO_MARKER_GENES_DESCRIPTION,
  TABLE_HEADER_SPECIFICITY,
  SPECIFICITY_TOOLTIP_CONTENT_FIRST_HALF,
  SPECIFICITY_TOOLTIP_CONTENT_SECOND_HALF,
  DIFFERENTIAL_EXPRESSION_LINK_TEXT,
} from "./constants";
import { useConnect } from "./connect";
import {
  DivTable,
  DivTableCell,
  DivTableCellPadded,
  DivTableHead,
  DivTableRow,
} from "../../common/styles";
import Description from "src/views/CellGuide/components/CellGuideCard/components/Description";
import { StyledQuestionMarkIcon } from "src/common/style";

function CellInfoSideBar({
  cellInfoCellType,
  tissueInfo,
  generateGeneInfo,
}: CellInfoBarProps): JSX.Element | null {
  const {
    isLoading,
    data,
    handleCopyGenes,
    handleDisplayGenes,
    handleFmgHoverEnd,
    handleMarkerScoreHoverEnd,
    handleSpecificityHoverEnd,
    setHoverStartTime,
    differentialExpressionUrl,
  } = useConnect({
    cellInfoCellType,
  });

  if (isLoading || !data) return null;

  if (!cellInfoCellType) return null;

  const numMarkerGenes = Object.keys(data.marker_genes).length;
  const shouldShowEmptyState =
    numMarkerGenes === 0 ||
    cellInfoCellType.cellType.total_count < 25 ||
    TISSUES_WITHOUT_MARKER_GENES.includes(tissueInfo.name);

  return (
    <>
      <TissueName>{tissueInfo.name}</TissueName>
      <Description
        cellTypeId={cellInfoCellType.cellType.id}
        cellTypeName={cellInfoCellType.cellType.name}
        skinnyMode={true}
        inSideBar
      />
      <StyledLink
        href={`${ROUTES.CELL_GUIDE}/${cellInfoCellType.cellType.id}`}
        onClick={() =>
          track(EVENTS.WMG_OPEN_IN_CG_CLICKED, {
            cell_type: cellInfoCellType.cellType.id,
          })
        }
        target="_blank"
        rel="noreferrer noopener"
      >
        {MARKER_SCORE_CELLGUIDE_LINK_TEXT}
        <Icon sdsIcon="ChevronRight" sdsType="static" sdsSize="xs" />
      </StyledLink>
      <StyledLink
        href={differentialExpressionUrl}
        onClick={() =>
          track(EVENTS.WMG_OPEN_IN_DE_CLICKED, {
            cell_type: cellInfoCellType.cellType.id,
            tissue: cellInfoCellType.tissueID,
          })
        }
        target="_blank"
        rel="noreferrer noopener"
      >
        {DIFFERENTIAL_EXPRESSION_LINK_TEXT}
        <Icon sdsIcon="ChevronRight" sdsType="static" sdsSize="xs" />
      </StyledLink>
      <ButtonContainer>
        <ButtonWrapper>
          <StyledMarkerGeneHeader>{MARKER_GENE_LABEL}</StyledMarkerGeneHeader>
          <Tooltip
            sdsStyle="dark"
            placement="bottom"
            width="default"
            className="fmg-tooltip-icon"
            arrow
            onOpen={() => setHoverStartTime(Date.now())}
            onClose={handleFmgHoverEnd}
            title={
              <StyledTooltip>
                <TooltipContent>{MARKER_GENES_TOOLTIP_CONTENT}</TooltipContent>
                <TooltipLink
                  href={ROUTES.FMG_DOCS}
                  rel="noopener"
                  target="_blank"
                  onClick={() => {
                    track(EVENTS.WMG_FMG_QUESTION_BUTTON_HOVER, {
                      label: MARKER_GENE_LABEL,
                    });
                    track(EVENTS.WMG_FMG_DOCUMENTATION_CLICKED, {
                      label: MARKER_GENE_LABEL,
                    });
                  }}
                >
                  {MARKER_SCORE_TOOLTIP_LINK_TEXT}
                </TooltipLink>
              </StyledTooltip>
            }
          >
            <TooltipButton
              sdsStyle="minimal"
              sdsType="secondary"
              isAllCaps={false}
            >
              <StyledQuestionMarkIcon />
            </TooltipButton>
          </Tooltip>
          <BetaChip label="Beta" size="small" />
        </ButtonWrapper>
        <AddToDotplotButton
          data-testid="add-to-dotplot-fmg-button"
          startIcon={<Icon sdsIcon="Plus" sdsSize="s" sdsType="button" />}
          onClick={handleDisplayGenes}
          sdsStyle="minimal"
          sdsType="primary"
          isAllCaps={false}
          disabled={shouldShowEmptyState}
        >
          {MARKER_SCORE_DOTPLOT_BUTTON_TEXT}
        </AddToDotplotButton>
      </ButtonContainer>
      {shouldShowEmptyState ? (
        (track(EVENTS.WMG_FMG_NO_MARKER_GENES, {
          combination: `${cellInfoCellType.cellType.id}, ${tissueInfo.id}`,
        }),
        (
          <NoMarkerGenesContainer data-testid="no-marker-genes-warning">
            <NoMarkerGenesHeader>{NO_MARKER_GENES_HEADER}</NoMarkerGenesHeader>
            <NoMarkerGenesDescription data-testid="no-marker-genes-description">
              {TISSUES_WITHOUT_MARKER_GENES.includes(tissueInfo.name)
                ? NO_MARKER_GENES_FOR_BLOOD_DESCRIPTION
                : cellInfoCellType.cellType.total_count < 25
                ? TOO_FEW_CELLS_NO_MARKER_GENES_DESCRIPTION
                : NO_MARKER_GENES_DESCRIPTION}
            </NoMarkerGenesDescription>
          </NoMarkerGenesContainer>
        ))
      ) : (
        <DivTable>
          <DivTableHead>
            <DivTableCell>
              {TABLE_HEADER_GENE}
              <CopyGenesButton
                onClick={handleCopyGenes}
                sdsType="primary"
                sdsStyle="minimal"
                isAllCaps={false}
                startIcon={<Icon sdsIcon="Copy" sdsSize="s" sdsType="button" />}
              >
                Copy
              </CopyGenesButton>
            </DivTableCell>
            <DivTableCell align data-testid="marker-genes-table-header-score">
              {EFFECT_SIZE}
              <Tooltip
                sdsStyle="dark"
                placement="bottom"
                width="default"
                className="fmg-tooltip-icon"
                arrow
                onOpen={() => setHoverStartTime(Date.now())}
                onClose={handleMarkerScoreHoverEnd}
                title={
                  <StyledTooltip>
                    <TooltipContent>
                      {MARKER_SCORE_TOOLTIP_CONTENT}
                      <br />
                      <br />
                      {MARKER_SCORE_LOW_CONTENT}
                      <br />
                      {MARKER_SCORE_MEDIUM_CONTENT}
                      <br />
                      {MARKER_SCORE_HIGH_CONTENT}
                      <br />
                    </TooltipContent>
                    <TooltipLink
                      href={ROUTES.FMG_DOCS}
                      rel="noopener"
                      target="_blank"
                      onClick={() => {
                        track(EVENTS.WMG_FMG_QUESTION_BUTTON_HOVER, {
                          label: MARKER_SCORE_LABEL,
                        });
                        track(EVENTS.WMG_FMG_DOCUMENTATION_CLICKED, {
                          label: MARKER_SCORE_LABEL,
                        });
                      }}
                    >
                      {MARKER_SCORE_TOOLTIP_LINK_TEXT}
                    </TooltipLink>
                  </StyledTooltip>
                }
              >
                <TooltipButton
                  sdsStyle="minimal"
                  sdsType="secondary"
                  isAllCaps={false}
                  data-testid={MARKER_SCORE_TOOLTIP_TEST_ID}
                >
                  <StyledQuestionMarkIcon />
                </TooltipButton>
              </Tooltip>
            </DivTableCell>
            <DivTableCell align data-testid="marker-genes-table-specificity">
              {TABLE_HEADER_SPECIFICITY}
              <Tooltip
                sdsStyle="dark"
                placement="bottom"
                width="default"
                className="fmg-tooltip-icon"
                arrow
                onOpen={() => setHoverStartTime(Date.now())}
                onClose={handleSpecificityHoverEnd}
                title={
                  <StyledTooltip>
                    <TooltipContent>
                      {SPECIFICITY_TOOLTIP_CONTENT_FIRST_HALF} {tissueInfo.name}{" "}
                      {SPECIFICITY_TOOLTIP_CONTENT_SECOND_HALF}
                    </TooltipContent>
                    <TooltipLink
                      href={ROUTES.FMG_DOCS}
                      rel="noopener"
                      target="_blank"
                      onClick={() => {
                        track(EVENTS.WMG_FMG_QUESTION_BUTTON_HOVER, {
                          label: MARKER_SCORE_LABEL,
                        });
                        track(EVENTS.WMG_FMG_DOCUMENTATION_CLICKED, {
                          label: MARKER_SCORE_LABEL,
                        });
                      }}
                    >
                      {MARKER_SCORE_TOOLTIP_LINK_TEXT}
                    </TooltipLink>
                  </StyledTooltip>
                }
              >
                <TooltipButton
                  sdsStyle="minimal"
                  sdsType="secondary"
                  isAllCaps={false}
                  data-testid={SPECIFICITY_TOOLTIP_TEST_ID}
                >
                  <StyledQuestionMarkIcon />
                </TooltipButton>
              </Tooltip>
            </DivTableCell>
          </DivTableHead>
          {Object.entries(data.marker_genes).map(([symbol, metadata]) => (
            <DivTableRow key={symbol}>
              <DivTableCell>
                {symbol}
                <InfoButtonWrapper
                  data-testid="gene-info-button-cell-info"
                  onClick={() => {
                    generateGeneInfo(symbol);

                    track(EVENTS.WMG_FMG_GENE_INFO, {
                      gene: symbol,
                    });
                  }}
                >
                  <InfoSVG id="gene-info-button-fmg" />
                </InfoButtonWrapper>
              </DivTableCell>
              <DivTableCellPadded data-testid="marker-scores-fmg" align>
                {metadata.marker_score.toFixed(2)}
              </DivTableCellPadded>
              <DivTableCellPadded data-testid="specificity-fmg" align>
                {metadata.specificity.toFixed(2)}
              </DivTableCellPadded>
            </DivTableRow>
          ))}
        </DivTable>
      )}
    </>
  );
}

export default React.memo(CellInfoSideBar);
