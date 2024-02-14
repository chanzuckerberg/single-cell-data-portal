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
        isAllCaps={false}
        onClick={onCancelReorder}
        sdsStyle="minimal"
        sdsType="secondary"
      >
        Cancel
      </MinimalButton>
      <SquareButton
        color="success"
        onClick={onSaveReorder}
        sdsStyle="square"
        sdsType="primary"
      >
        Save Order
      </SquareButton>
    </>
  );
}
