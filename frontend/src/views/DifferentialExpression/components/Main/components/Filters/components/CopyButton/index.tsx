import { QueryGroup } from "src/views/DifferentialExpression/common/store/reducer";
import { StyledButtonIcon } from "./style";
import { useConnect } from "./connect";
interface Props {
  queryGroupKey: keyof QueryGroup;
  testId?: string;
}
function CopyButton({ queryGroupKey, testId }: Props): JSX.Element {
  const { handleClick, disabled } = useConnect({ queryGroupKey });
  return (
    <StyledButtonIcon
      onClick={handleClick}
      sdsIcon="copy"
      sdsSize="small"
      disabled={disabled}
      data-testid={testId}
    />
  );
}
export default CopyButton;
