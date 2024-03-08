export const getFormattedExplorerUrl = ({
  selectedOrganism,
  tissueId,
  cellTypeId,
}: {
  selectedOrganism?: string;
  tissueId: string;
  cellTypeId?: string;
}): { explorerUrl: string; formattedSelectedOrganism: string } => {
  // construct explorer URL
  const formattedSelectedOrganism =
    selectedOrganism?.toLowerCase()?.replace(/ /g, "_") ?? "";
  const formattedTissueId = tissueId.replace(/:/g, "_");
  const formattedCellTypeId = cellTypeId?.replace(/:/g, "_");
  let explorerUrl = "";

  if (formattedTissueId && formattedCellTypeId && formattedSelectedOrganism) {
    explorerUrl = `https://cellxgene.cziscience.com/e/cellguide-cxgs/tissues/${formattedSelectedOrganism}/${formattedTissueId}__${formattedCellTypeId}.cxg/`;
  } else if (formattedCellTypeId && formattedSelectedOrganism) {
    explorerUrl = `https://cellxgene.cziscience.com/e/cellguide-cxgs/${formattedSelectedOrganism}/${formattedCellTypeId}.cxg/`;
  } else if (formattedTissueId && formattedSelectedOrganism) {
    explorerUrl = `https://cellxgene.cziscience.com/e/cellguide-cxgs/tissues/${formattedSelectedOrganism}/${formattedTissueId}__CL_0000000.cxg/`;
  }
  return {
    explorerUrl,
    formattedSelectedOrganism,
  };
};
