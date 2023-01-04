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
import React, { useCallback, useContext } from "react";
import { useMarkerGenes } from "src/common/queries/wheresMyGene";
import { BetaChip } from "src/components/Header/style";
import { DispatchContext, State } from "../../common/store";
import { addSelectedGenes } from "../../common/store/actions";
import {
  ButtonContainer,
  CopyGenesButton,
  GeneCellHeader,
  GeneHeaderWrapper,
  TissueName,
  TooltipButton,
} from "./style";
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
  }, [data]);

  const handleDisplayGenes = useCallback(() => {
    if (!data || !dispatch) return;
    const genes = Object.keys(data.marker_genes);
    dispatch(addSelectedGenes(genes));
  }, [data, dispatch]);

  if (isLoading || !data) return null;

  if (!cellInfoCellType) return null;
  return (
    <div>
      <TissueName>{tissueName}</TissueName>
      <ButtonContainer>
        <div>
          <Tooltip
            sdsStyle="dark"
            placement="bottom"
            width="default"
            arrow={true}
            title={
              <>
                <div>Marker genes are highly and uniquely expressed in the cell type relative to all other cell types.</div>
                <br/>
                <div><a href={``} style={{color: "inherit"}}>Click to read more about the identification method.</a></div>
              </>
            }
          >
          <TooltipButton
            endIcon={<Icon sdsIcon="infoCircle" sdsSize="s" sdsType="button" />}
            onClick={handleCopyGenes}
            sdsStyle="minimal"
            sdsType="secondary"
            isAllCaps={false}
            style={{ fontWeight: "500" }}
          >
            Marker Genes
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
            P-value
          </CellHeader>
          <CellHeader hideSortIcon horizontalAlign="right">
            Effect Size
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
                primaryText={gene[1].p_value.toPrecision(4)}
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
