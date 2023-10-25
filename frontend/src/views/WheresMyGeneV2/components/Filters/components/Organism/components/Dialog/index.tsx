import { Dialog as RawDialog } from "@czi-sds/components";
import {
  ActionButton,
  StyledDialogAction,
  StyledDialogContent,
  StyledDialogPaper,
  StyledDialogTitle,
  Title,
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
      <StyledDialogTitle
        title={(<Title>Change organism?</Title>) as unknown as string}
        data-testid="dialog-title"
      />
      <StyledDialogContent data-testid="dialog-content">
        This will reset the dot plot and remove all selected genes and filters.
      </StyledDialogContent>
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
