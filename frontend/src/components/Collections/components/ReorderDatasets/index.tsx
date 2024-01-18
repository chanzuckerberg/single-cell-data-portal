import {
  MinimalButton,
  SquareButton,
} from "src/components/Collections/components/ReorderDatasets/style";
import { OnReorderFn, REORDER_MODE } from "src/common/hooks/useReorderMode";

interface Props {
  onReorder: OnReorderFn;
}

export default function ReorderDatasets({ onReorder }: Props): JSX.Element {
  return (
    <>
      <MinimalButton
        isAllCaps={false}
        onClick={() => onReorder(REORDER_MODE.INACTIVE)}
        sdsStyle="minimal"
        sdsType="secondary"
      >
        Cancel
      </MinimalButton>
      <SquareButton
        color="success"
        onClick={() => onReorder(REORDER_MODE.INACTIVE)} // TODO(cc) implement save function.
        sdsStyle="square"
        sdsType="primary"
      >
        Save Order
      </SquareButton>
    </>
  );
}
