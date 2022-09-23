import { Row } from "react-table";
import {
  Categories,
  TableCountSummary,
} from "src/components/common/Filter/common/entities";
import { DATASET_MAX_CELL_COUNT } from "src/components/common/Grid/common/constants";

/**
 * Returns table count summary.
 * @param filteredRows - Filtered result set.
 * @param originalRows - Original result set before filtering.
 * @returns Table count summary with filtered rows count and total rows count.
 */
export function buildTableCountSummary<T extends Categories>(
  filteredRows: Row<T>[],
  originalRows: Row<T>[]
): TableCountSummary {
  let rowCount = 0;
  let totalCount = 0;
  if (filteredRows?.length && originalRows?.length) {
    rowCount = filteredRows.length;
    totalCount = originalRows.length;
  }
  return { row: rowCount, total: totalCount };
}

/**
 * Returns true if cell count is over the max allowable cell count.
 * @param cellCount
 * @returns true if cell count is over the max allowable cell count
 */
export function checkIsOverMaxCellCount(cellCount: number | null): boolean {
  return (cellCount || 0) > DATASET_MAX_CELL_COUNT;
}
