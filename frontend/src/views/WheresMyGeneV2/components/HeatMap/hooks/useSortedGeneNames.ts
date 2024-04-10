import { agnes } from "ml-hclust";
import { useMemo } from "react";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import {
  CellType,
  CellTypeGeneExpressionSummaryData,
  GeneExpressionSummary,
  SORT_BY,
  Tissue,
} from "src/views/WheresMyGeneV2/common/types";

export const TISSUE_CELL_TYPE_DIVIDER = "~";

interface Props {
  tissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap: Map<
    string,
    Map<string, Map<string, CellTypeGeneExpressionSummaryData>>
  >;
  selectedCellTypes: { [tissue: Tissue]: CellType[] };
  geneSortBy: SORT_BY;
  genes: string[];
}

export function useSortedGeneNames({
  tissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap,
  selectedCellTypes,
  geneSortBy,
  genes,
}: Props): string[] {
  const isSortByUserEntered = geneSortBy === SORT_BY.USER_ENTERED;

  const geneNamesJSON = useMemo(() => {
    return JSON.stringify(genes);
  }, [genes]);

  const columns: string[] = useMemo(() => {
    if (isSortByUserEntered) return [];

    const result = new Set<string>();

    for (const [tissueName, cellTypes] of Object.entries(selectedCellTypes)) {
      for (const cellType of cellTypes) {
        result.add(
          `${tissueName}${TISSUE_CELL_TYPE_DIVIDER}${cellType.viewId}`
        );
      }
    }

    return Array.from(result);
  }, [selectedCellTypes, isSortByUserEntered]);

  const orderedGeneNamesJSON = useMemo(() => {
    if (isSortByUserEntered) return geneNamesJSON;

    const matrix = genes.map((geneName) => {
      return columns.map((column) => {
        const [tissueName, cellTypeId] = column.split(TISSUE_CELL_TYPE_DIVIDER);

        const cellTypeIdToGeneNameToCellTypeGeneExpressionSummaryData =
          tissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap.get(
            tissueName
          );

        if (!cellTypeIdToGeneNameToCellTypeGeneExpressionSummaryData) {
          return 0;
        }

        const geneNameToCellTypeGeneExpressionSummaryData =
          cellTypeIdToGeneNameToCellTypeGeneExpressionSummaryData.get(
            cellTypeId
          );

        if (!geneNameToCellTypeGeneExpressionSummaryData) return 0;

        const geneExpressionSummaryData =
          geneNameToCellTypeGeneExpressionSummaryData.get(geneName);

        if (!geneExpressionSummaryData) return 0;

        const { meanExpression = 0, percentage = 0 } =
          geneExpressionSummaryData;

        return meanExpression * percentage;
      });
    });

    const tree = agnes(matrix);

    const orderedGeneNames = tree?.indices().map((index) => genes[index]);

    return JSON.stringify(orderedGeneNames || EMPTY_ARRAY);
  }, [
    columns,
    genes,
    tissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap,
    isSortByUserEntered,
    geneNamesJSON,
  ]);

  // (thuang): Return cached value if JSON values are the same
  return useMemo(() => {
    return geneSortBy === SORT_BY.USER_ENTERED
      ? JSON.parse(geneNamesJSON)
      : JSON.parse(orderedGeneNamesJSON);
  }, [geneSortBy, geneNamesJSON, orderedGeneNamesJSON]);
}

export function useTissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap(selectedGeneExpressionSummariesByTissueName: {
  [tissueName: string]: GeneExpressionSummary[];
}): Map<string, Map<string, Map<string, CellTypeGeneExpressionSummaryData>>> {
  return useMemo(() => {
    return getTissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap(
      selectedGeneExpressionSummariesByTissueName
    );
  }, [selectedGeneExpressionSummariesByTissueName]);
}

function getTissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap(selectedGeneExpressionSummariesByTissueName: {
  [tissueName: string]: GeneExpressionSummary[];
}) {
  const result = new Map<
    Tissue,
    Map<string, Map<string, CellTypeGeneExpressionSummaryData>>
  >();

  for (const [tissueName, geneExpressionSummaries] of Object.entries(
    selectedGeneExpressionSummariesByTissueName
  )) {
    const cellTypeIdToGeneNameToCellTypeGeneExpressionSummaryData =
      result.get(tissueName) ||
      new Map<string, Map<string, CellTypeGeneExpressionSummaryData>>();

    for (const geneExpressionSummary of geneExpressionSummaries) {
      const { name: geneName, cellTypeGeneExpressionSummaries = [] } =
        geneExpressionSummary || {};

      if (!geneName || !cellTypeGeneExpressionSummaries.length) continue;

      for (const cellTypeGeneExpressionSummary of cellTypeGeneExpressionSummaries) {
        const geneNameToCellTypeGeneExpressionSummaryData =
          cellTypeIdToGeneNameToCellTypeGeneExpressionSummaryData.get(
            cellTypeGeneExpressionSummary.viewId
          ) || new Map<string, CellTypeGeneExpressionSummaryData>();

        geneNameToCellTypeGeneExpressionSummaryData.set(
          geneName,
          cellTypeGeneExpressionSummary
        );

        cellTypeIdToGeneNameToCellTypeGeneExpressionSummaryData.set(
          cellTypeGeneExpressionSummary.viewId,
          geneNameToCellTypeGeneExpressionSummaryData
        );
      }
    }

    result.set(
      tissueName,
      cellTypeIdToGeneNameToCellTypeGeneExpressionSummaryData
    );
  }

  return result;
}
