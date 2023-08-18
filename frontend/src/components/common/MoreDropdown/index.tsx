import { IButtonProps, IPopoverProps, Popover } from "@blueprintjs/core";
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
        data-testid="collection-more-button"
        sdsIcon="dotsHorizontal"
        sdsSize="small"
        sdsType="tertiary"
      />
    </Popover>
  );
};

export default MoreDropdown;
