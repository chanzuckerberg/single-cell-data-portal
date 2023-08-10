import { Button, Icon, Tooltip } from "@czi-sds/components";
import React, { useCallback, useContext, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { ROUTES } from "src/common/constants/routes";
import { OntologyTerm, useMarkerGenes } from "src/common/queries/wheresMyGene";
import { BetaChip } from "src/components/Header/style";
import { DispatchContext, State } from "../../common/store";
import { addSelectedGenes } from "../../common/store/actions";
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
  TissueName,
  TooltipButton,
} from "./style";
import { Link } from "../../../../components/GeneInfoSideBar/style";
import questionMarkIcon from "src/common/images/question-mark-icon.svg";
import { StyledImage } from "../HeatMap/components/YAxisChart/style";
import InfoSVG from "../HeatMap/components/YAxisChart/icons/info-sign-icon.svg";
import { RightSidebarProperties } from "../../../../components/common/RightSideBar";
import { InfoButtonWrapper } from "src/components/common/Filter/common/style";

const MARKER_GENE_LABEL = "marker genes";
const MARKER_SCORE_LABEL = "marker score";

export interface CellInfoBarProps extends RightSidebarProperties {
  cellInfoCellType: Exclude<State["cellInfoCellType"], null>;
  tissueInfo: OntologyTerm;
  generateGeneInfo: (gene: string) => void;
}

function CellInfoSideBar({
  cellInfoCellType,
  tissueInfo,
  generateGeneInfo,
}: CellInfoBarProps): JSX.Element | null {
  const urlParams = new URLSearchParams(window.location.search);
  let testType: "binomtest" | undefined = undefined;

  if (urlParams.get("test") === "binomtest") {
    testType = "binomtest";
  }
  const { isLoading, data } = useMarkerGenes({
    cellTypeID: cellInfoCellType.cellType.id,
    organismID: cellInfoCellType.organismID,
    test: testType,
    tissueID: cellInfoCellType.tissueID,
  });

  const dispatch = useContext(DispatchContext);

  const handleCopyGenes = useCallback(() => {
    if (!data) return;
    const genes = Object.keys(data.marker_genes);
    navigator.clipboard.writeText(genes.join(", "));
    track(EVENTS.WMG_FMG_COPY_GENES_CLICKED);
  }, [data]);

  const handleDisplayGenes = useCallback(() => {
    if (!data || !dispatch) return;
    const genes = Object.keys(data.marker_genes);
    dispatch(addSelectedGenes(genes));
    track(EVENTS.WMG_FMG_ADD_GENES_CLICKED);
  }, [data, dispatch]);

  const [hoverStartTime, setHoverStartTime] = useState(0);

  const useHandleHoverEnd = (event: EVENTS, payload = {}) => {
    return useCallback(() => {
      if (Date.now() - hoverStartTime > 2 * 1000) {
        track(event, payload);
      }
    }, [hoverStartTime]);
  };

  const handleFmgHoverEnd = useHandleHoverEnd(
    EVENTS.WMG_FMG_QUESTION_BUTTON_HOVER,
    { label: MARKER_GENE_LABEL }
  );
  const handleMarkerScoreHoverEnd = useHandleHoverEnd(
    EVENTS.WMG_FMG_QUESTION_BUTTON_HOVER,
    { label: MARKER_SCORE_LABEL }
  );

  if (isLoading || !data) return null;

  const numMarkerGenes = Object.keys(data.marker_genes).length;

  if (!cellInfoCellType) return null;
  return (
    <div>
      <TissueName>{tissueInfo.name}</TissueName>
      <Link
        href={`${ROUTES.CELL_GUIDE}/${cellInfoCellType.cellType.id}`}
        target="_blank"
        rel="noreferrer noopener"
      >
        Open in CellGuide
      </Link>
      <ButtonContainer>
        <div>
          <StyledMarkerGeneHeader>Marker Genes</StyledMarkerGeneHeader>
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
                <div>
                  Marker genes are highly and uniquely expressed in the cell
                  type relative to all other cell types.
                </div>
                <br />
                <div>
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
                    Click to read more about the identification method.
                  </a>
                </div>
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
        </div>
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
          Add to Dot Plot
        </Button>
      </ButtonContainer>

      {!numMarkerGenes ? (
        (track(EVENTS.WMG_FMG_NO_MARKER_GENES, {
          combination: `${cellInfoCellType.cellType.id}, ${tissueInfo.id}`,
        }),
        (
          <NoMarkerGenesContainer data-testid="no-marker-genes-warning">
            <NoMarkerGenesHeader>No Marker Genes</NoMarkerGenesHeader>
            <NoMarkerGenesDescription>
              No reliable marker genes for this cell type.
            </NoMarkerGenesDescription>
          </NoMarkerGenesContainer>
        ))
      ) : (
        <StyledHTMLTable condensed bordered={false}>
          <thead>
            <tr>
              <td>Gene </td>
              <td>
                Marker Score
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
                      <div>
                        Marker Score indicates the strength and specificity of a
                        gene as a marker. It is the 5th percentile of the effect
                        sizes when comparing the expressions in a cell type of
                        interest to each other cell type in the tissue.
                      </div>
                      <br />
                      <div>
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
                          Click to read more about the identification method.
                        </a>
                      </div>
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
    </div>
  );
}

export default React.memo(CellInfoSideBar);
