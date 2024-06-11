import { StyledButton } from "src/components/Collections/components/ReorderDatasets/style";
import { ReorderAction } from "src/views/Collection/hooks/useReorder/useReorder";
import { Button } from "src/components/common/Button";

interface Props {
  reorderAction: ReorderAction;
}

export default function ReorderDatasets({ reorderAction }: Props): JSX.Element {
  const { onCancelReorder, onSaveReorder } = reorderAction;
  return (
    <>
      <StyledButton
        data-testid="datasets-reorder-cancel"
        isAllCaps={false}
        onClick={onCancelReorder}
        sdsStyle="minimal"
        sdsType="secondary"
      >
        Cancel
      </StyledButton>
      <Button
        data-testid="datasets-reorder-save"
        color="success"
        onClick={onSaveReorder}
        sdsStyle="square"
        sdsType="primary"
      >
        Save Order
      </Button>
    </>
  );
}
