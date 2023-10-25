import { Dialog as RawDialog, DialogContent } from "@czi-sds/components";
import {
  ActionButton,
  StyledDialogAction,
  StyledDialogPaper,
  StyledDialogTitle,
} from "./style";
import { Button } from "src/components/common/Button";

interface Props {
  handleCancel: () => void;
  handleConfirm: () => void;
  isOpen: boolean;
}

export default function Dialog({ handleCancel, handleConfirm, isOpen }: Props) {
  return (
    <RawDialog
      PaperComponent={StyledDialogPaper}
      sdsSize="xs"
      canClickOutsideClose={false}
      open={isOpen}
    >
      <StyledDialogTitle title="Change organism?" data-testid="dialog-title" />
      <DialogContent data-testid="dialog-content">
        This will reset the dot plot and remove all selected genes and filters.
      </DialogContent>
      <StyledDialogAction data-testid="dialog-actions" buttonPosition="right">
        <>
          <Button
            isAllCaps={false}
            sdsStyle="minimal"
            sdsType="secondary"
            onClick={handleCancel}
          >
            Cancel
          </Button>
          <ActionButton
            color="error"
            sdsStyle="square"
            sdsType="primary"
            onClick={handleConfirm}
          >
            Confirm
          </ActionButton>
        </>
      </StyledDialogAction>
    </RawDialog>
  );
}
