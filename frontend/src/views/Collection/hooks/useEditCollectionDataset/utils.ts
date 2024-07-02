import {
  EditDataset,
  UseEditCollectionDataset,
} from "src/views/Collection/hooks/useEditCollectionDataset/types";
import { CollectionNotification } from "src/views/Collection/hooks/useNotification/types";
import { EDIT_DATASET_ERROR_NOTIFICATION } from "src/views/Collection/hooks/useEditCollectionDataset/constants";

/**
 * Returns edit dataset related actions related to dataset editing.
 * @param editCollectionDataset - Edit collection dataset.
 * @param collectionNotification - Collection notification.
 * @return edit dataset related actions.
 */
export function getEditDataset(
  editCollectionDataset: UseEditCollectionDataset,
  collectionNotification: CollectionNotification
): EditDataset {
  const { editDatasetAction, formMethod } = editCollectionDataset;
  const { onEditDataset } = editDatasetAction;
  const { onNotify } = collectionNotification;
  return {
    onEditDataset,
    onError: () => onNotify(EDIT_DATASET_ERROR_NOTIFICATION),
    formMethod,
  };
}
