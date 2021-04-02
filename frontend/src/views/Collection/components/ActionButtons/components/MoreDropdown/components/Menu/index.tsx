import {
  IMenuItemProps,
  Intent,
  Menu as RawMenu,
  MenuItem,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import React from "react";
import DeleteCollection from "src/components/Collections/components/DeleteCollection";

const DeleteButton = (props: IMenuItemProps) => {
  return (
    <MenuItem
      {...props}
      shouldDismissPopover={false}
      icon={IconNames.TRASH}
      intent={Intent.DANGER}
      text="Delete Collection"
    />
  );
};

interface Props {
  id?: string;
}

const Menu = ({ id = "" }: Props) => {
  return (
    <RawMenu>
      <DeleteCollection id={id} Button={DeleteButton} />
    </RawMenu>
  );
};

export default Menu;
