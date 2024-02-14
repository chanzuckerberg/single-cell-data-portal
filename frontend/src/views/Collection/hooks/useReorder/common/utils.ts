import { UseReorder } from "src/views/Collection/hooks/useReorder/useReorder";
import { Reorder } from "src/views/Collection/hooks/useReorder/common/entities";
import {
  Dataset,
  PROCESSING_STATUS,
  UPLOAD_STATUS,
  VALIDATION_STATUS,
} from "src/common/entities";

/**
 * Returns the datasets IDs, associated with the given collection.
 * @param datasets - Datasets.
 * @returns collection dataset IDs.
 */
export function getDatasetIds(datasets: Dataset[]): string[] {
  return datasets.map((dataset) => dataset.id);
}

/**
 * Returns reorder related values (i.e. disabled) and actions related to dataset reordering.
 * @param reorder - Reorder.
 * @param datasets - Datasets.
 * @return reorder related values and actions.
 */
export function getReorder(reorder: UseReorder, datasets: Dataset[]): Reorder {
  const { reorderAction } = reorder;
  const { onStartReorder } = reorderAction;
  return {
    ...reorder,
    disabled: !isCollectionDatasetsReorderable(datasets),
    startReorder: () => onStartReorder(getDatasetIds(datasets)),
  };
}

/**
 * Returns true if the collection datasets are reorder-able (i.e. there are more than one dataset, and all datasets are uploaded and valid).
 * @param datasets - Datasets.
 * @returns true if the collection datasets are reorder-able.
 */
export function isCollectionDatasetsReorderable(datasets: Dataset[]): boolean {
  return (
    datasets.length > 1 && datasets.every(isCollectionDatasetUploadedAndValid)
  );
}

/**
 * Returns true if the dataset is uploaded and valid, with a processing status of success.
 * @param dataset - Dataset.
 * @returns true if the dataset is uploaded and valid.
 */
function isCollectionDatasetUploadedAndValid(dataset: Dataset): boolean {
  return (
    dataset.processing_status.processing_status === PROCESSING_STATUS.SUCCESS &&
    dataset.processing_status.upload_status === UPLOAD_STATUS.UPLOADED &&
    dataset.processing_status.validation_status === VALIDATION_STATUS.VALID
  );
}
