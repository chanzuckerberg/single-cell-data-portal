import styled from "@emotion/styled";
import { Button, Icon } from "czifui";
import React, { useCallback, useContext } from "react";
import { useMarkerGenes } from "src/common/queries/wheresMyGene";
import { BetaChip } from "src/components/Header/style";
import { DispatchContext, State } from "../../common/store";
import { addSelectedGenes } from "../../common/store/actions";
export interface CellInfoBarProps {
  cellInfoCellTypes: State["cellInfoCellTypes"];
}
function CellInfoBar({
  cellInfoCellTypes,
}: CellInfoBarProps): JSX.Element | null {
  const { isLoading, data } = useMarkerGenes({
    cellTypeID: cellInfoCellTypes[0].cellType.id,
    tissueID: cellInfoCellTypes[0].tissueID,
    organismID: cellInfoCellTypes[0].organismID,
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

  const ButtonContainer = styled.div`
    display: flex;
    flex-direction: row;
  `;

  const cellInfoCellType = cellInfoCellTypes[0];
  return (
    <div>
      <h3>{cellInfoCellType.cellType.name}</h3>
      <BetaChip label="Beta" size="small" />
      <ButtonContainer>
        <Button
          startIcon={<Icon sdsIcon="infoCircle" sdsSize="s" sdsType="button" />}
          onClick={handleCopyGenes}
          sdsStyle="minimal"
          sdsType="info"
          isAllCaps={false}
        >
          Marker Genes
        </Button>
        <Button
          startIcon={<Icon sdsIcon="copy" sdsSize="s" sdsType="button" />}
          onClick={handleCopyGenes}
          sdsStyle="minimal"
          sdsType="primary"
          isAllCaps={false}
        >
          Copy Genes
        </Button>
        <Button
          startIcon={<Icon sdsIcon="plus" sdsSize="s" sdsType="button" />}
          onClick={handleDisplayGenes}
          sdsStyle="minimal"
          sdsType="primary"
          isAllCaps={false}
        >
          Add to Dot Plot
        </Button>
      </ButtonContainer>
      <br />
      <table>
        <thead>
          <tr>
            <td>Gene</td>
            <td>P-value</td>
            <td>Effect Size</td>
          </tr>
        </thead>
        {Object.entries(data.marker_genes).map((gene) => (
          <tr key={gene[0]}>
            <td>{gene[0]}</td>
            <td>{gene[1].p_value}</td>
            <td>{gene[1].effect_size}</td>
          </tr>
        ))}
      </table>
    </div>
  );
}

export default CellInfoBar;
