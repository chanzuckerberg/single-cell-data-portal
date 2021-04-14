import {
  IMenuItemProps,
  Intent,
  Menu as RawMenu,
  MenuItem,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import React from "react";
import DeleteCollection from "src/components/Collections/components/DeleteCollection";
import CreateCollection from "src/components/CreateCollectionModal";
import styled from "styled-components";

const DeleteButton = (props: Partial<IMenuItemProps>) => {
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

const EditButton = (props: Partial<IMenuItemProps>) => {
  return (
    <MenuItem
      {...props}
      shouldDismissPopover={false}
      icon={IconNames.EDIT}
      intent={Intent.NONE}
      text="Edit Details"
    />
  );
};

interface Props {
  id?: string;
}

const StyledMenu = styled(RawMenu)`
  border-radius: 3px;
  padding: 8px;
  & > li:last-child {
    margin-bottom: 0;
  }
`;

const Menu = ({ id = "" }: Props) => {
  return (
    <StyledMenu>
      <DeleteCollection id={id} Button={DeleteButton} />
      <EditButton />
      <CreateCollection id={id} />
    </StyledMenu>
  );
};

export default Menu;
