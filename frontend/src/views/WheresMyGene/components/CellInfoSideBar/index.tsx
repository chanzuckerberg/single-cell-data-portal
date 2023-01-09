import {
  Button,
  CellBasic,
  CellHeader,
  Icon,
  Table,
  TableHeader,
  TableRow,
  Tooltip,
} from "czifui";
import React, { useCallback, useContext, useEffect, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { ROUTES } from "src/common/constants/routes";
import { useMarkerGenes } from "src/common/queries/wheresMyGene";
import { BetaChip } from "src/components/Header/style";
import { DispatchContext, State } from "../../common/store";
import { addSelectedGenes } from "../../common/store/actions";
import {
  ButtonContainer,
  CopyGenesButton,
  GeneCellHeader,
  GeneHeaderWrapper,
  StyledIconImage,
  StyledMarkerGeneHeader,
  StyledTooltip,
  TissueName,
  TooltipButton,
} from "./style";
import questionMarkIcon from "src/common/images/question-mark-icon.svg";
export interface CellInfoBarProps {
  cellInfoCellType: Exclude<State["cellInfoCellType"], null>;
  tissueName: string;
}

function CellInfoSideBar({
  cellInfoCellType,
  tissueName,
}: CellInfoBarProps): JSX.Element | null {
  const urlParams = new URLSearchParams(window.location.search);
  let testType: "ttest" | undefined = undefined;

  if (urlParams.get("test") === "ttest") {
    testType = "ttest";
  }
  const { isLoading, data } = useMarkerGenes({
    cellTypeID: cellInfoCellType.cellType.id,
    tissueID: cellInfoCellType.tissueID,
    organismID: cellInfoCellType.organismID,
    test: testType,
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

  const [hoverTime, setHoverTime] = useState({ start: -1, end: -1 });

  useEffect(() => {
    const { start, end } = hoverTime;
    if (start !== -1 && end !== -1 && end - start > 2 * 1000) {
      track(EVENTS.WMG_FMG_QUESTION_BUTTON_HOVER);
      setHoverTime({ start: -1, end: -1 });
    }
  }, [hoverTime]);

  if (isLoading || !data) return null;

  if (!cellInfoCellType) return null;
  return (
    <div>
      <TissueName>{tissueName}</TissueName>
      <ButtonContainer>
        <div>
          <StyledMarkerGeneHeader>Marker Genes</StyledMarkerGeneHeader>
          <Tooltip
            sdsStyle="dark"
            placement="bottom"
            width="default"
            className="fmg-tooltip-icon"
            arrow={true}
            title={
              <StyledTooltip>
                <div>Marker genes are highly and uniquely expressed in the cell type relative to all other cell types.</div>
                <br/>
                <div>
                  <a 
                    href={ROUTES.FMG_DOCS} 
                    rel="noopener" 
                    target="_blank" 
                    onClick={() => track(EVENTS.WMG_FMG_DOCUMENTATION_CLICKED)}
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
              style={{ fontWeight: "500" }}
              onHoverStart={() => setHoverTime({ start: Date.now(), end: -1 })}
              onHoverEnd={() =>
                setHoverTime({ start: hoverTime.start, end: Date.now() })
              }
            >
              <StyledIconImage src={questionMarkIcon} />
            </TooltipButton>
          </Tooltip>
          <BetaChip label="Beta" size="small" />
        </div>
        <Button
          startIcon={<Icon sdsIcon="plus" sdsSize="s" sdsType="button" />}
          onClick={handleDisplayGenes}
          sdsStyle="minimal"
          sdsType="primary"
          isAllCaps={false}
          style={{ fontWeight: "500" }}
        >
          Add to Dot Plot
        </Button>
      </ButtonContainer>
      <Table>
        <TableHeader>
          <GeneCellHeader hideSortIcon>
            <GeneHeaderWrapper>
              Gene{" "}
              <CopyGenesButton
                onClick={handleCopyGenes}
                sdsType="primary"
                sdsStyle="minimal"
                isAllCaps={false}
                startIcon={<Icon sdsIcon="copy" sdsSize="s" sdsType="button" />}
              >
                Copy
              </CopyGenesButton>
            </GeneHeaderWrapper>
          </GeneCellHeader>
          <CellHeader hideSortIcon horizontalAlign="right">
            Marker Score
          </CellHeader>
        </TableHeader>
        <tbody>
          {Object.entries(data.marker_genes).map((gene) => (
            <TableRow key={gene[0]}>
              <CellBasic
                shouldShowTooltipOnHover={false}
                primaryText={gene[0]}
              />
              <CellBasic
                shouldShowTooltipOnHover={false}
                horizontalAlign="right"
                primaryText={gene[1].effect_size.toPrecision(4)}
              />
            </TableRow>
          ))}
        </tbody>
      </Table>
    </div>
  );
}

export default CellInfoSideBar;
