import { Props } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/types";
import React, { Fragment } from "react";
import { Dialog as StyledDialog } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/style";
import { DialogActions, DialogContent, DialogTitle } from "@czi-sds/components";
import { Button as StyledButton } from "src/components/common/Button";
import {
  CANCEL_BUTTON_PROPS,
  DIALOG_PROPS,
  SAVE_BUTTON_PROPS,
} from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/constants";
import { useDialog } from "src/views/Collection/hooks/useDialog";
import { useEditCollectionDataset } from "src/views/Collection/hooks/useEditCollectionDataset/useEditCollectionDataset";
import EditDatasetFields from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/components/EditDatasetFields";
import { useForm } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/hooks/useForm";
import { Dataset } from "src/common/entities";
import { FieldValues } from "src/views/Collection/hooks/useEditCollectionDataset/common/entities";

export default function EditDataset({
  Button,
  collectionId,
  dataset,
}: Props): JSX.Element {
  const { onClose, onOpen, open } = useDialog();
  const {
    editDatasetAction: { onEditDataset },
  } = useEditCollectionDataset();
  const { clearErrors, errors, handleSubmit } = useForm();
  const fieldValues = mapDatasetToFieldValues(dataset);
  const { id: datasetId } = dataset;
  return (
    <Fragment>
      <Button onClick={onOpen} />
      <StyledDialog
        {...DIALOG_PROPS}
        onClose={() => onClose(clearErrors)}
        onSubmit={handleSubmit(
          onEditDataset,
          { collectionId, datasetId },
          fieldValues
        )}
        open={open}
      >
        <DialogTitle title="Edit Dataset" />
        <DialogContent>
          <EditDatasetFields errors={errors} fieldValues={fieldValues} />
        </DialogContent>
        <DialogActions>
          <StyledButton
            {...CANCEL_BUTTON_PROPS}
            onClick={() => onClose(clearErrors)}
          >
            Cancel
          </StyledButton>
          <StyledButton {...SAVE_BUTTON_PROPS}>Save</StyledButton>
        </DialogActions>
      </StyledDialog>
    </Fragment>
  );
}

/**
 * Maps dataset to form field values.
 * @param dataset - Dataset.
 * @returns dataset mapped to form field values.
 */
function mapDatasetToFieldValues(dataset: Dataset): FieldValues {
  return {
    title: dataset.name,
  };
}
