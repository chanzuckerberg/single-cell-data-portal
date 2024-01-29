import {
  MinimalButton,
  SquareButton,
} from "src/components/Collections/components/ReorderDatasets/style";
import {
  OnSetReorderModeFn,
  REORDER_MODE,
} from "src/common/hooks/useReorderMode";

interface Props {
  onSetReorderMode: OnSetReorderModeFn;
}

export default function ReorderDatasets({
  onSetReorderMode,
}: Props): JSX.Element {
  return (
    <>
      <MinimalButton
        isAllCaps={false}
        onClick={() => onSetReorderMode(REORDER_MODE.INACTIVE)}
        sdsStyle="minimal"
        sdsType="secondary"
      >
        Cancel
      </MinimalButton>
      <SquareButton
        color="success"
        onClick={() => onSetReorderMode(REORDER_MODE.INACTIVE)} // TODO(cc) implement save function.
        sdsStyle="square"
        sdsType="primary"
      >
        Save Order
      </SquareButton>
    </>
  );
}
