import { AnchorButton } from "@blueprintjs/core";
import { ReactElement } from "react";
import { ActionButton as StyledActionButton } from "src/components/common/Grid/components/ActionButton/style";

interface Props {
  iconSvg: ReactElement;
  isAnchorButton?: boolean;
}

export default function ActionButton({
  iconSvg,
  isAnchorButton = false,
  ...props /* Spread props to allow for data-test-id and button specific attributes e.g. "href", "target", or "disabled". */
}: Props): JSX.Element {
  return (
    <StyledActionButton
      as={isAnchorButton ? AnchorButton : undefined}
      icon={iconSvg}
      minimal
      {...props}
    />
  );
}
