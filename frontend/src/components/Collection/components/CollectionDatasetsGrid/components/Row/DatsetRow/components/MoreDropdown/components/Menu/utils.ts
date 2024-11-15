import { Dataset } from "src/common/entities";

/**
 * Determines whether the delete dataset menu option should be available.
 * A dataset may be deleted if the collection is private, or the dataset has not been previously
 * published; where the published_at property is used to determine whether the dataset has been previously published.
 * @param dataset - Dataset.
 * @param revisionsEnabled - Dataset revisions enabled.
 * @returns true if the delete dataset menu option should be available.
 */
export function isDeleteDatasetAvailable(
  dataset: Dataset,
  revisionsEnabled: boolean
): boolean {
  return !revisionsEnabled || !dataset.published_at;
}
