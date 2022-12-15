import { Button, Icon } from "czifui";
import React, { useCallback, useContext } from "react";
import { useMarkerGenes } from "src/common/queries/wheresMyGene";
import { BetaChip } from "src/components/Header/style";
import { DispatchContext, State } from "../../common/store";
import { addSelectedGenes } from "../../common/store/actions";
import {
  ButtonContainer,
  CopyGenesButton,
  StyledHTMLTable,
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
      <StyledHTMLTable condensed bordered={false}>
        <thead>
          <tr>
            <td>
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
            </td>
            <td>P-value</td>
            <td>Effect Size</td>
          </tr>
        </thead>
        <tbody>
          {Object.entries(data.marker_genes).map((gene) => (
            <tr key={gene[0]}>
              <td>{gene[0]}</td>
              <td>{gene[1].p_value.toPrecision(4)}</td>
              <td>{gene[1].effect_size.toPrecision(4)}</td>
            </tr>
          ))}
        </tbody>
      </StyledHTMLTable>
    </div>
  );
}

export default CellInfoSideBar;
