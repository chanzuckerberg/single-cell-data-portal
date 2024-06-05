import { Dialog as RawDialog } from "@czi-sds/components";
import {
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
            onClick={handleCancel}
            sdsStyle="minimal"
            sdsType="secondary"
          >
            Cancel
          </Button>
          <Button
            color="error"
            onClick={handleConfirm}
            sdsStyle="square"
            sdsType="destructive"
          >
            Confirm
          </Button>
        </>
      </StyledDialogAction>
    </RawDialog>
  );
}
