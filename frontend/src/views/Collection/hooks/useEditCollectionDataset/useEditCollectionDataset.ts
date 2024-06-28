import {
  FieldValues,
  UseEditCollectionDataset,
} from "src/views/Collection/hooks/useEditCollectionDataset/common/entities";
import { Dataset } from "src/common/entities";
import { useEditDataset } from "src/common/queries/collections";

/**
 * Edit functionality for collection dataset.
 */
export function useEditCollectionDataset(): UseEditCollectionDataset {
  const editDatasetMutation = useEditDataset();

  const onEditDataset = async (
    dataset: Dataset,
    fieldValues: FieldValues
  ): Promise<void> => {
    // Send order to BE.
    const { collection_id: collectionId, id: datasetId } = dataset;
    const payload = JSON.stringify(fieldValues);
    await editDatasetMutation.mutateAsync(
      {
        collectionId,
        datasetId,
        payload,
      },
      {
        onSuccess: () => {
          // TODO(cc) close dialog.
        },
        onError: () => {
          // TODO(cc) error banner.
        },
      }
    );
  };

  return {
    editDatasetAction: { onEditDataset },
  };
}
