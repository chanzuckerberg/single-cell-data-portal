import { Cluster, agnes } from "ml-hclust";
import { useMemo } from "react";
import { EMPTY_ARRAY, noop } from "src/common/constants/utils";
import { useEffectDebugger } from "src/common/hooks/useEffectDebugger";
import { CellTypeRow } from "src/common/queries/wheresMyGene";
import {
  CellTypeGeneExpressionSummaryData,
  SORT_BY,
  Tissue,
} from "src/views/WheresMyGeneV2/common/types";

interface Props {
  tissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap: Map<
    string,
    Map<string, Map<string, CellTypeGeneExpressionSummaryData>>
  >;
  selectedCellTypes: { [tissue: Tissue]: CellTypeRow[] };
  genes: string[];
  cellTypeSortBy: SORT_BY;
}

export function useSortedCellTypesByTissueName({
  tissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap,
  selectedCellTypes,
  genes,
  cellTypeSortBy,
}: Props): { [tissue: Tissue]: CellTypeRow[] } {
  const conditionalDependencies =
    cellTypeSortBy === SORT_BY.CELL_ONTOLOGY
      ? EMPTY_ARRAY
      : [
          genes,
          tissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap,
        ];

  return useMemo(() => {
    console.log("recomputing sorted cell types");
    const isSortByCellOntology = cellTypeSortBy === SORT_BY.CELL_ONTOLOGY;

    if (isSortByCellOntology) {
      return selectedCellTypes;
    }

    const sortedCellTypesByTissueName: { [tissueName: string]: CellTypeRow[] } =
      {};

    for (const [tissueName, cellTypes] of Object.entries(selectedCellTypes)) {
      const cellTypeIdToGeneNameToCellTypeGeneExpressionSummaryData =
        tissueNameToCellTypeIdToGeneNameToCellTypeGeneExpressionSummaryDataMap.get(
          tissueName
        );

      if (!cellTypeIdToGeneNameToCellTypeGeneExpressionSummaryData) continue;

      const { aggregatedOnly, optionOnly } = separateCellTypes(cellTypes);

      const cellTypeOptionsByCellTypeId =
        buildCellTypeOptionsByCellTypeId(optionOnly);

      const tree = hierarchicalClustering({
        aggregatedOnly,
        cellTypeIdToGeneNameToCellTypeGeneExpressionSummaryData,
        genes,
      });

      const orderedCellTypes = [];

      for (const index of tree.indices()) {
        const cellType = aggregatedOnly[index];

        const cellTypeOptions = cellTypeOptionsByCellTypeId[cellType.id] || [];

        // (thuang): Reverse the order, so `unknown` shows up last
        if (cellType.id.includes("UBERON"))
          orderedCellTypes.unshift(cellType, ...cellTypeOptions.reverse());
        else orderedCellTypes.push(cellType, ...cellTypeOptions.reverse());
      }

      /**
       * (thuang): Reverse the order so the first item shows up at the top of
       * the heatmap
       */
      sortedCellTypesByTissueName[tissueName] = orderedCellTypes.reverse();
    }
    return sortedCellTypesByTissueName;
    // eslint-disable-next-line react-hooks/exhaustive-deps -- conditional dependencies are dependent on cellTypeSortBy
  }, [selectedCellTypes, cellTypeSortBy, ...conditionalDependencies]);
}

function hierarchicalClustering({
  aggregatedOnly,
  cellTypeIdToGeneNameToCellTypeGeneExpressionSummaryData,
  genes,
}: {
  aggregatedOnly: CellTypeRow[];
  cellTypeIdToGeneNameToCellTypeGeneExpressionSummaryData: Map<
    string,
    Map<string, CellTypeGeneExpressionSummaryData>
  >;
  genes: string[];
}): Cluster {
  const matrix = aggregatedOnly.map((cellType) => {
    const geneNameToCellTypeGeneExpressionSummaryData =
      cellTypeIdToGeneNameToCellTypeGeneExpressionSummaryData.get(
        cellType.viewId
      );

    if (!geneNameToCellTypeGeneExpressionSummaryData) {
      return genes.map(() => 0);
    }

    return genes.map((geneName) => {
      const cellTypeGeneExpressionSummaryData =
        geneNameToCellTypeGeneExpressionSummaryData.get(geneName);

      if (!cellTypeGeneExpressionSummaryData) return 0;

      const { meanExpression, percentage } = cellTypeGeneExpressionSummaryData;

      if (percentage === undefined) return 0;

      return meanExpression * percentage;
    });
  });

  return agnes(matrix);
}

function separateCellTypes(cellTypes: CellTypeRow[]): {
  aggregatedOnly: CellTypeRow[];
  optionOnly: CellTypeRow[];
} {
  const aggregatedOnly: CellTypeRow[] = [];
  const optionOnly: CellTypeRow[] = [];

  cellTypes.forEach((cellType) => {
    if (cellType.isAggregated) {
      aggregatedOnly.push(cellType);
    } else {
      optionOnly.push(cellType);
    }
  });

  return { aggregatedOnly, optionOnly };
}

function buildCellTypeOptionsByCellTypeId(optionOnly: CellTypeRow[]) {
  const cellTypeOptionsByCellTypeId: {
    [cellTypeId: CellTypeRow["id"]]: CellTypeRow[];
  } = {};

  for (const cellType of optionOnly) {
    const cellTypeOptions = cellTypeOptionsByCellTypeId[cellType.id] || [];

    cellTypeOptions.push(cellType);

    cellTypeOptionsByCellTypeId[cellType.id] = cellTypeOptions;
  }

  return cellTypeOptionsByCellTypeId;
}
