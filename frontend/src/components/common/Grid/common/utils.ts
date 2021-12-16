import { DATASET_MAX_CELL_COUNT } from "src/components/common/Grid/common/constants";

/**
 * Returns true if cell count is over the max allowable cell count.
 * @param cellCount
 * @returns true if cell count is over the max allowable cell count
 */
export function checkIsOverMaxCellCount(cellCount: number | null): boolean {
  return (cellCount || 0) > DATASET_MAX_CELL_COUNT;
}
