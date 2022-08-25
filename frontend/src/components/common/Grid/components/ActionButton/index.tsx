import { ReactElement } from "react";
import {
  StyledAnchorButton,
  StyledButton,
} from "src/components/common/Grid/components/ActionButton/style";

interface Props {
  iconSvg: ReactElement;
  isAnchorButton?: boolean;
}

export default function ActionButton({
  iconSvg,
  isAnchorButton = false,
  ...props /* Spread props to allow for data-test-id and button specific attributes e.g. "href", "target", or "disabled". */
}: Props): JSX.Element {
  const Component = isAnchorButton ? StyledAnchorButton : StyledButton;

  return <Component icon={iconSvg} minimal {...props} />;
}
