import { Dataset } from "src/common/entities";

/**
 * TODO(cc) remove when dataset order is available.
 * @param {Dataset[]} datasets
 * @return {Dataset[]}
 */
export function sortByCellCountDescending(datasets: Dataset[]): Dataset[] {
  return (
    datasets?.sort((a, b) => (b.cell_count ?? 0) - (a.cell_count ?? 0)) || []
  );
}

/**
 * Returns sorted datasets based on the given dataset order.
 * @param datasets - Datasets.
 * @param datasetIDs - Dataset IDs, ordered.
 * @returns sorted datasets.
 */
export function sortByDatasetIDOrder(
  datasets: Dataset[],
  datasetIDs: string[]
): Dataset[] {
  return datasetIDs
    .map((datasetID) => datasets.find(({ id }) => id === datasetID))
    .filter(Boolean) as Dataset[];
}
