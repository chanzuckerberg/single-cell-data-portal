import { Row } from "react-table";
import {
  Categories,
  TableCountSummary,
} from "src/components/common/Filter/common/entities";
import { DATASET_MAX_CELL_COUNT } from "src/components/common/Grid/common/constants";

/**
 * Basic compare function, returning a sort return value.
 * See React Table https://github.com/TanStack/table/blob/beccddcab001434f3bb11843b3fda72f8b000cc2/packages/table-core/src/sortingFns.ts#L53
 * @param val0 - First value.
 * @param val1 - Second value.
 * @returns sort return value (0 | 1 | -1).
 */
function basicSort<TValue>(val0: TValue, val1: TValue): number {
  return val0 === val1 ? 0 : val0 > val1 ? 1 : -1;
}

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

/**
 * Returns sort return value from the compare function for row values of type "array".
 * Uses a compare function based off React Table "basic" compare function.
 * See React Table https://github.com/TanStack/table/blob/beccddcab001434f3bb11843b3fda72f8b000cc2/packages/table-core/src/sortingFns.ts.
 * @param rowA - First row to sort.
 * @param rowB - Second row to sort.
 * @param columnId - Sorted column identifier.
 * @returns sort return value from the compare function (0 | 1 | -1).
 */
export function arraySortingFn<T extends object>(
  rowA: Row<T>,
  rowB: Row<T>,
  columnId: string
): number {
  const rowAValue = rowA.values[columnId];
  const rowBValue = rowB.values[columnId];
  const rowALength = rowAValue?.length || 0;
  const rowBLength = rowBValue?.length || 0;
  // Start with a comparison between the first value of each array and continue until a non-zero result is found.
  for (let i = 0; i < Math.min(rowALength, rowBLength); i++) {
    const result = basicSort(
      toString(rowAValue[i]).toLowerCase(),
      toString(rowBValue[i]).toLowerCase()
    );
    if (result !== 0) {
      return result; // Return the first non-zero result.
    }
  }
  // If each array value is equal, then compare the length of the arrays.
  // The array with the least values is sorted first.
  return Math.sign(rowALength - rowBLength);
}

/**
 * Returns the given value as a string.
 * See React Table "toString" function
 * https://github.com/TanStack/table/blob/beccddcab001434f3bb11843b3fda72f8b000cc2/packages/table-core/src/sortingFns.ts#L57
 * @param tValue - Cell value.
 * @returns the value as a string.
 */
function toString<TValue>(tValue: TValue): string {
  if (typeof tValue === "number") {
    if (isNaN(tValue) || tValue === Infinity || tValue === -Infinity) {
      return "";
    }
    return String(tValue);
  }
  if (typeof tValue === "string") {
    return tValue;
  }
  return "";
}
