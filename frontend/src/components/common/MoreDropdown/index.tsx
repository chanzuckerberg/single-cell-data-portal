import {
  Button,
  IButtonProps,
  IPopoverProps,
  Popover,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import React from "react";

interface Props {
  popoverProps?: IPopoverProps;
  buttonProps?: IButtonProps;
}

const MoreDropdown = ({ popoverProps = {}, buttonProps = {} }: Props) => {
  return (
    <Popover {...popoverProps}>
      <Button {...buttonProps} minimal icon={IconNames.MORE} />
    </Popover>
  );
};

export default MoreDropdown;
