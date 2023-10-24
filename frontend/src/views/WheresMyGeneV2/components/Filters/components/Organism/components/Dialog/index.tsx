import {
  Dialog as RawDialog,
  DialogContent,
  DialogActions,
} from "@czi-sds/components";
import { ActionButton, CancelButton, StyledDialogTitle } from "./style";

interface Props {
  handleCancel: () => void;
  handleConfirm: () => void;
  isOpen: boolean;
}

export default function Dialog({ handleCancel, handleConfirm, isOpen }: Props) {
  return (
    <RawDialog sdsSize="xs" canClickOutsideClose={false} open={isOpen}>
      <StyledDialogTitle title="Change organism?" data-testid="dialog-title" />
      <DialogContent data-testid="dialog-content">
        This will reset the dot plot and remove all selected genes and filters.
      </DialogContent>
      <DialogActions data-testid="dialog-actions" buttonPosition="right">
        <>
          <CancelButton
            isAllCaps={false}
            sdsStyle="minimal"
            sdsType="secondary"
            onClick={handleCancel}
          >
            Cancel
          </CancelButton>
          <ActionButton
            color="error"
            sdsStyle="square"
            sdsType="primary"
            onClick={handleConfirm}
          >
            Confirm
          </ActionButton>
        </>
      </DialogActions>
    </RawDialog>
  );
}
