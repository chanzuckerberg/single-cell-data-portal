import {
  Button,
  IButtonProps,
  IPopoverProps,
  Popover,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";

interface Props {
  popoverProps?: IPopoverProps;
  buttonProps?: IButtonProps;
}

const MoreDropdown = ({
  popoverProps = {},
  buttonProps = {},
}: Props): JSX.Element => {
  return (
    <Popover {...popoverProps}>
      <Button
        {...buttonProps}
        minimal
        icon={IconNames.MORE}
        data-test-id="collection-more-button"
      />
    </Popover>
  );
};

export default MoreDropdown;
