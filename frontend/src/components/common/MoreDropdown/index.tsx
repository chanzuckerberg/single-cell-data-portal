import { IButtonProps, IPopoverProps, Popover } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { MoreButton } from "src/components/common/MoreDropdown/style";

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
      <MoreButton
        {...buttonProps}
        minimal
        icon={IconNames.MORE}
        data-testid="collection-more-button"
      />
    </Popover>
  );
};

export default MoreDropdown;
