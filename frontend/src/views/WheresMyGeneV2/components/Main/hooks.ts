import { useMemo } from "react";

import {
  CellTypeByTissueName,
  GeneExpressionSummariesByTissueName,
  getOntologyTermIdFromCellTypeViewId,
} from "src/common/queries/wheresMyGene";

import { State } from "src/views/WheresMyGeneV2/common/store";

export function useScaledMeanExpression({
  cellTypesByTissueName,
  geneExpressionSummariesByTissueName,
  selectedGenes,
  expandedTissueIds,
  filteredCellTypes,
  filteredTissueIds,
}: {
  cellTypesByTissueName: CellTypeByTissueName;
  geneExpressionSummariesByTissueName: GeneExpressionSummariesByTissueName;
  selectedGenes: State["selectedGenes"];
  expandedTissueIds: State["expandedTissueIds"];
  filteredCellTypes: State["filteredCellTypes"];
  filteredTissueIds: State["selectedFilters"]["tissues"];
}) {
  // TODO(thuang): Fix this complexity
  // eslint-disable-next-line sonarjs/cognitive-complexity
  return useMemo(() => {
    let min = Infinity;
    let max = -Infinity;
    for (const [tissueName, tissueSelectedCellTypes] of Object.entries(
      cellTypesByTissueName
    )) {
      const tissueId = tissueSelectedCellTypes.at(-1)?.id ?? "";
      const tissueSelectedCellTypeIds = tissueSelectedCellTypes.map(
        (cellType) => cellType.viewId
      );
      const tissueSelectedCellTypeNames = tissueSelectedCellTypes.map(
        (cellType) => cellType.cellTypeName
      );
      // get object of cell type id to cell type name
      const cellTypeNameById: { [id: string]: string } = {};
      tissueSelectedCellTypes.forEach((cellType) => {
        cellTypeNameById[cellType.id] = cellType.cellTypeName;
      });
      const tissueGeneExpressionSummaries =
        geneExpressionSummariesByTissueName[tissueName];

      if (!tissueGeneExpressionSummaries) {
        continue;
      }

      for (const selectedGeneName of selectedGenes) {
        const geneExpressionSummary =
          tissueGeneExpressionSummaries[selectedGeneName];

        if (geneExpressionSummary) {
          const { cellTypeGeneExpressionSummaries } = geneExpressionSummary;

          for (const cellTypeGeneExpressionSummary of cellTypeGeneExpressionSummaries) {
            // get term before $
            const cellTypeGeneExpressionSummaryId =
              getOntologyTermIdFromCellTypeViewId(
                cellTypeGeneExpressionSummary.viewId
              );
            const isCellType =
              !cellTypeGeneExpressionSummaryId.startsWith("UBERON");

            /*
            This is to dynamically set the min and max based on the data that is visible to users.
            Skip conditions:
            - If the cell type is not part of an expanded tissue
            - If the tissue does not contain any of the filtered cell types and at least one cell type is filtered
            - If the cell type is not included in the filtered cell types and at least one cell type is filtered
            - If the tissue id is not included in the filtered tissue ids and at least one tissue is filtered
            - If the cell type gene expression summary view id is not included in the available view ids
            The above conditions capture all possible scenarios where the data is not visible to users.
            */
            if (
              // If the cell type is not part of an expanded tissue
              (isCellType && !expandedTissueIds.includes(tissueId)) ||
              // If the tissue does not contain any of the filtered cell types and at least one cell type is filtered
              (!isCellType &&
                filteredCellTypes.length > 0 &&
                !tissueSelectedCellTypeNames.filter((value) =>
                  filteredCellTypes.includes(value)
                ).length) ||
              // If the cell type is not included in the filtered cell types and at least one cell type is filtered
              (isCellType &&
                filteredCellTypes.length > 0 &&
                !filteredCellTypes.includes(
                  cellTypeNameById?.[cellTypeGeneExpressionSummaryId]
                )) ||
              // If the tissue id is not included in the filtered tissue ids and at least one tissue is filtered
              (!isCellType &&
                filteredTissueIds.length > 0 &&
                !filteredTissueIds.includes(tissueId)) ||
              // If the cell type gene expression summary view id is not included in the available view ids
              !tissueSelectedCellTypeIds.includes(
                cellTypeGeneExpressionSummary.viewId
              )
            ) {
              continue;
            }

            const { meanExpression } = cellTypeGeneExpressionSummary;

            min = Math.min(min, meanExpression);
            max = Math.max(max, meanExpression);
          }
        }
      }
    }
    return {
      scaledMeanExpressionMax: max,
      scaledMeanExpressionMin: min,
    };
  }, [
    geneExpressionSummariesByTissueName,
    cellTypesByTissueName,
    selectedGenes,
    filteredTissueIds,
    expandedTissueIds,
    filteredCellTypes,
  ]);
}
