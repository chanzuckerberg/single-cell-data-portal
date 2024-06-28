import { Popover, PopoverProps } from "@blueprintjs/core";
import { Button, ButtonProps } from "@czi-sds/components";

interface Props {
  popoverProps?: PopoverProps;
  buttonProps?: Partial<ButtonProps>;
}

const MoreDropdown = ({
  popoverProps = {},
  buttonProps = {},
}: Props): JSX.Element => {
  return (
    <Popover {...popoverProps}>
      <Button
        {...buttonProps}
        data-testid="collection-more-button"
        icon="DotsHorizontal"
        sdsSize="small"
        sdsStyle="icon"
        sdsType="secondary"
      />
    </Popover>
  );
};

export default MoreDropdown;
