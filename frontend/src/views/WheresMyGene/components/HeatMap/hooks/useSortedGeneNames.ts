import { agnes } from "ml-hclust";
import { useMemo } from "react";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { State } from "src/views/WheresMyGene/common/store";
import {
  CellType,
  CellTypeGeneExpressionSummaryData,
  Genes,
  SORT_BY,
  Tissue,
} from "src/views/WheresMyGene/common/types";
import { SelectedGeneExpressionSummariesByTissueName } from "..";

const TISSUE_CELL_TYPE_DIVIDER = "~";

interface Props {
  tissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap: Map<
    string,
    Map<string, Map<string, CellTypeGeneExpressionSummaryData>>
  >;
  selectedCellTypes: { [tissue: Tissue]: CellType[] };
  geneSortBy: SORT_BY;
  genes: State["selectedGenes"];
}

export function useSortedGeneNames({
  tissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap,
  selectedCellTypes,
  geneSortBy,
  genes,
}: Props): {[groupName: string]: string[]} {
  const isSortByUserEntered = geneSortBy === SORT_BY.USER_ENTERED;

  const geneNamesJSON = useMemo(() => {
    return JSON.stringify(Object.fromEntries(genes));
  }, [genes]);

  const columns: string[] = useMemo(() => {
    if (isSortByUserEntered) return [];

    const result = new Set<string>();

    for (const [tissueName, cellTypes] of Object.entries(selectedCellTypes)) {
      for (const cellType of cellTypes) {
        result.add(`${tissueName}${TISSUE_CELL_TYPE_DIVIDER}${cellType.id}`);
      }
    }

    return Array.from(result);
  }, [selectedCellTypes, isSortByUserEntered]);

  const orderedGeneNamesJSON = useMemo(() => {
    if (isSortByUserEntered) return geneNamesJSON;

    const orderedGeneNamesByGroupName: {[groupName: string]: Genes} = {};
    for (const [groupName, geneList] of genes) {
      const matrix = geneList.map((geneName) => {
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

      const orderedGeneNames = tree?.indices().map((index) => geneList[index]);  
      orderedGeneNamesByGroupName[groupName] = orderedGeneNames; 
    }




    return JSON.stringify(orderedGeneNamesByGroupName || EMPTY_ARRAY);
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

export function useTissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap(
  selectedGeneExpressionSummariesByTissueName: SelectedGeneExpressionSummariesByTissueName): Map<string, Map<string, Map<string, CellTypeGeneExpressionSummaryData>>> {
  return useMemo(() => {
    return getTissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap(
      selectedGeneExpressionSummariesByTissueName
    );
  }, [selectedGeneExpressionSummariesByTissueName]);
}

function getTissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap(
  selectedGeneExpressionSummariesByTissueName: SelectedGeneExpressionSummariesByTissueName) {
  const result = new Map<
    Tissue,
    Map<string, Map<string, CellTypeGeneExpressionSummaryData>>
  >();

  for (const [_, geneExpressionSummariesByTissue] of Object.entries(
    selectedGeneExpressionSummariesByTissueName
  )) {
    for (const [tissueName, geneExpressionSummaries] of Object.entries(
      geneExpressionSummariesByTissue
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
              cellTypeGeneExpressionSummary.id
            ) || new Map<string, CellTypeGeneExpressionSummaryData>();

          geneNameToCellTypeGeneExpressionSummaryData.set(
            geneName,
            cellTypeGeneExpressionSummary
          );

          cellTypeIdToGeneNameToCellTypeGeneExpressionSummaryData.set(
            cellTypeGeneExpressionSummary.id,
            geneNameToCellTypeGeneExpressionSummaryData
          );
        }
      }

      result.set(
        tissueName,
        cellTypeIdToGeneNameToCellTypeGeneExpressionSummaryData
      );      
    }
  }

  return result;
}
