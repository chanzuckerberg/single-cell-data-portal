/**
 * Helper functions for binding dataset metadata.
 */

// App dependencies
import { Dataset } from "../../common/types/entities";

/**
 * Filter out datasets that have cell counts greater than the given cell count limit.
 * @param datasets - Array of datasets to filter.
 * @param cellCountLimit - Maximum number of cells a dataset can have in order to be included for display.
 * @returns Updated metadata with datasets with more than the given cell count limit filtered out.
 */
export function removeLargeDatasets(
  datasets: Dataset[],
  cellCountLimit: number,
): Dataset[] {
  return datasets.filter((dataset) => dataset.cell_count <= cellCountLimit);
}
