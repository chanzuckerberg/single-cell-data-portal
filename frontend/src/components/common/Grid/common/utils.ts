import { DatasetDeployment } from "src/common/entities";
import { DATASET_MAX_CELL_COUNT } from "src/components/common/Grid/common/constants";

/**
 * Returns true if cell count is over the max allowable cell count.
 * @param cellCount
 * @returns true if cell count is over the max allowable cell count
 */
export function checkIsOverMaxCellCount(cellCount: number | null): boolean {
  return (cellCount || 0) > DATASET_MAX_CELL_COUNT;
}

/**
 * Returns true if dataset deployments are valid.
 * TODO(cc) copy and modified from frontend/src/components/Collections/components/Grid/components/Row/DatasetRow/utils.ts
 * @param dataset_deployments
 */
export function hasCXGFile(dataset_deployments: DatasetDeployment[]): boolean {
  const deployments = dataset_deployments;

  if (!deployments || !deployments.length) return false;

  return true;
}
