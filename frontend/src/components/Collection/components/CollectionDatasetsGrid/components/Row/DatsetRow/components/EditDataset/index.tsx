import { Props } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/types";
import React, { Fragment, useCallback } from "react";
import { Dialog as StyledDialog } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/style";
import { DialogActions, DialogContent, DialogTitle } from "@czi-sds/components";
import { Button as StyledButton } from "src/components/common/Button";
import {
  CANCEL_BUTTON_PROPS,
  DATASET_EDIT_CANCEL,
  DATASET_EDIT_SAVE,
  DIALOG_PROPS,
  SAVE_BUTTON_PROPS,
} from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/constants";
import { useDialog } from "src/views/Collection/hooks/useDialog";
import EditDatasetFields from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/components/EditDatasetFields";

export default function EditDataset({
  Button,
  collectionId,
  dataset,
  editDataset,
  menuProps,
}: Props): JSX.Element {
  const { onClose, onOpen, open } = useDialog();
  const { onEditDataset, formMethod } = editDataset;
  const {
    formState: { isDirty },
    handleSubmit,
    reset: resetForm,
  } = formMethod;
  const { id: datasetId } = dataset;

  const onAfterClose = useCallback(() => {
    resetForm();
    menuProps.onClose(); // Close the "more" menu after the dialog closes.
  }, [menuProps, resetForm]);

  return (
    <Fragment>
      <Button onClick={onOpen} />
      <StyledDialog
        {...DIALOG_PROPS}
        onClose={onClose}
        onSubmit={handleSubmit(
          onEditDataset,
          { collectionId, datasetId },
          { onError: onClose, onSuccess: onClose }
        )}
        onTransitionExited={onAfterClose}
        open={open}
      >
        <DialogTitle title="Edit Dataset" />
        <DialogContent>
          <EditDatasetFields dataset={dataset} formMethod={formMethod} />
        </DialogContent>
        <DialogActions>
          <StyledButton
            {...CANCEL_BUTTON_PROPS}
            data-testid={DATASET_EDIT_CANCEL}
            onClick={onClose}
          >
            Cancel
          </StyledButton>
          <StyledButton
            {...SAVE_BUTTON_PROPS}
            data-testid={DATASET_EDIT_SAVE}
            disabled={!isDirty}
          >
            Save
          </StyledButton>
        </DialogActions>
      </StyledDialog>
    </Fragment>
  );
}
