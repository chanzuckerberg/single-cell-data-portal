import { UseEditCollectionDataset } from "src/views/Collection/hooks/useEditCollectionDataset/types";
import { useEditDataset } from "src/common/queries/collections";
import {
  FieldValues,
  PathParams,
  SubmitOptions,
} from "src/views/Collection/hooks/useCollectionDatasetForm/types";
import { useCallback } from "react";
import { useCollectionDatasetForm } from "src/views/Collection/hooks/useCollectionDatasetForm/useCollectionDatasetForm";

/**
 * Edit functionality for collection dataset.
 */
export function useEditCollectionDataset(): UseEditCollectionDataset {
  const editDatasetFormMethod = useCollectionDatasetForm();
  const editDatasetMutation = useEditDataset();

  const onEditDataset = useCallback(
    async (
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
    },
    [editDatasetMutation]
  );

  return {
    editDatasetAction: { onEditDataset },
    formMethod: editDatasetFormMethod,
  };
}
