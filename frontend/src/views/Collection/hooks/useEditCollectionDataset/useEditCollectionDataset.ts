import {
  FieldValues,
  UseEditCollectionDataset,
} from "src/views/Collection/hooks/useEditCollectionDataset/common/entities";
import { useEditDataset } from "src/common/queries/collections";
import {
  PathParams,
  SubmitOptions,
} from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/hooks/types";

/**
 * Edit functionality for collection dataset.
 */
export function useEditCollectionDataset(): UseEditCollectionDataset {
  const editDatasetMutation = useEditDataset();

  const onEditDataset = async (
    pathParams: PathParams,
    fieldValues: FieldValues,
    submitOptions?: SubmitOptions
  ): Promise<void> => {
    // Send order to BE.
    const { collectionId, datasetId } = pathParams;
    const payload = JSON.stringify(fieldValues);
    await editDatasetMutation.mutateAsync(
      {
        collectionId,
        datasetId,
        payload,
      },
      {
        onSuccess: () => {
          submitOptions?.onSuccess?.();
        },
        onError: () => {
          submitOptions?.onError?.();
        },
      }
    );
  };

  return {
    editDatasetAction: { onEditDataset },
  };
}
