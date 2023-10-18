import { Button, Icon, Tooltip } from "@czi-sds/components";
import React from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { ROUTES } from "src/common/constants/routes";
import { BetaChip } from "src/components/Header/style";
import {
  ButtonContainer,
  CopyGenesButton,
  MarkerStrengthContainer,
  MarkerStrengthLabel,
  NoMarkerGenesContainer,
  NoMarkerGenesDescription,
  NoMarkerGenesHeader,
  StyledHTMLTable,
  StyledIconImage,
  StyledMarkerGeneHeader,
  StyledTooltip,
  TooltipContent,
  TissueName,
  TooltipButton,
} from "./style";
import { Link } from "../../../../components/GeneInfoSideBar/style";
import questionMarkIcon from "src/common/images/question-mark-icon.svg";
import { StyledImage } from "../HeatMap/components/YAxisChart/style";
import InfoSVG from "../HeatMap/components/YAxisChart/icons/info-sign-icon.svg";
import { InfoButtonWrapper } from "src/components/common/Filter/common/style";
import { CellInfoBarProps } from "./types";
import {
  MARKER_GENES_TOOLTIP_CONTENT,
  MARKER_GENE_LABEL,
  MARKER_SCORE_CELLGUIDE_LINK_TEXT,
  MARKER_SCORE_DOTPLOT_BUTTON_TEXT,
  MARKER_SCORE_LABEL,
  MARKER_SCORE_TOOLTIP_CONTENT,
  MARKER_SCORE_TOOLTIP_LINK_TEXT,
  NO_MARKER_GENES_DESCRIPTION,
  NO_MARKER_GENES_HEADER,
  TABLE_HEADER_GENE,
  TABLE_HEADER_SCORE,
} from "./constants";
import { useConnect } from "./connect";

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
    setHoverStartTime,
  } = useConnect({
    cellInfoCellType,
  });

  if (isLoading || !data) return null;

  const numMarkerGenes = Object.keys(data.marker_genes).length;

  if (!cellInfoCellType) return null;

  return (
    <>
      <TissueName>{tissueInfo.name}</TissueName>
      <Link
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
      </Link>
      <ButtonContainer>
        <>
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
                <>
                  <a
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
                  </a>
                </>
              </StyledTooltip>
            }
          >
            <TooltipButton
              sdsStyle="minimal"
              sdsType="secondary"
              isAllCaps={false}
            >
              <StyledIconImage src={questionMarkIcon} />
            </TooltipButton>
          </Tooltip>
          <BetaChip label="Beta" size="small" />
        </>
        <Button
          data-testid="add-to-dotplot-fmg-button"
          startIcon={<Icon sdsIcon="plus" sdsSize="s" sdsType="button" />}
          onClick={handleDisplayGenes}
          sdsStyle="minimal"
          sdsType="primary"
          isAllCaps={false}
          style={{ fontWeight: "500" }}
          disabled={!numMarkerGenes}
        >
          {MARKER_SCORE_DOTPLOT_BUTTON_TEXT}{" "}
        </Button>
      </ButtonContainer>
      {!numMarkerGenes ? (
        (track(EVENTS.WMG_FMG_NO_MARKER_GENES, {
          combination: `${cellInfoCellType.cellType.id}, ${tissueInfo.id}`,
        }),
        (
          <NoMarkerGenesContainer data-testid="no-marker-genes-warning">
            <NoMarkerGenesHeader>{NO_MARKER_GENES_HEADER}</NoMarkerGenesHeader>
            <NoMarkerGenesDescription>
              {NO_MARKER_GENES_DESCRIPTION}
            </NoMarkerGenesDescription>
          </NoMarkerGenesContainer>
        ))
      ) : (
        <StyledHTMLTable condensed bordered={false}>
          <thead>
            <tr>
              <td>{TABLE_HEADER_GENE} </td>
              <td>
                {TABLE_HEADER_SCORE}
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
                      </TooltipContent>
                      <>
                        <a
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
                        </a>
                      </>
                    </StyledTooltip>
                  }
                >
                  <TooltipButton
                    sdsStyle="minimal"
                    sdsType="secondary"
                    isAllCaps={false}
                  >
                    <StyledIconImage src={questionMarkIcon} />
                  </TooltipButton>
                </Tooltip>
              </td>
            </tr>
            <tr>
              <td>
                <CopyGenesButton
                  onClick={handleCopyGenes}
                  sdsType="primary"
                  sdsStyle="minimal"
                  isAllCaps={false}
                  startIcon={
                    <Icon sdsIcon="copy" sdsSize="s" sdsType="button" />
                  }
                >
                  Copy
                </CopyGenesButton>
              </td>
              <td>
                <MarkerStrengthContainer>
                  <MarkerStrengthLabel>{"Low: <1"}</MarkerStrengthLabel>
                  <MarkerStrengthLabel>{"Medium: 1-2"}</MarkerStrengthLabel>
                  <MarkerStrengthLabel>{"High: >2"}</MarkerStrengthLabel>
                </MarkerStrengthContainer>
              </td>
            </tr>
          </thead>
          <tbody>
            {Object.entries(data.marker_genes).map(([symbol, metadata]) => (
              <tr key={symbol}>
                <td>
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
                    <StyledImage
                      id="gene-info-button-fmg"
                      src={InfoSVG.src}
                      width="10"
                      height="10"
                      alt={`display gene info for ${symbol}`}
                    />
                  </InfoButtonWrapper>
                </td>
                <td data-testid="marker-scores-fmg">
                  {metadata.effect_size.toPrecision(4)}
                </td>
              </tr>
            ))}
          </tbody>
        </StyledHTMLTable>
      )}
    </>
  );
}

export default React.memo(CellInfoSideBar);
