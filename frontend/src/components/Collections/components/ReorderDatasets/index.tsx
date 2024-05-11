import {
  MinimalButton,
  SquareButton,
} from "src/components/Collections/components/ReorderDatasets/style";
import { ReorderAction } from "src/views/Collection/hooks/useReorder/useReorder";

interface Props {
  reorderAction: ReorderAction;
}

export default function ReorderDatasets({ reorderAction }: Props): JSX.Element {
  const { onCancelReorder, onSaveReorder } = reorderAction;
  return (
    <>
      <MinimalButton
        data-testid="datasets-reorder-cancel"
        isAllCaps={false}
        onClick={onCancelReorder}
        sdsType="secondary"
      >
        Cancel
      </MinimalButton>
      <SquareButton
        data-testid="datasets-reorder-save"
        color="success"
        onClick={onSaveReorder}
        sdsType="primary"
      >
        Save Order
      </SquareButton>
    </>
  );
}
