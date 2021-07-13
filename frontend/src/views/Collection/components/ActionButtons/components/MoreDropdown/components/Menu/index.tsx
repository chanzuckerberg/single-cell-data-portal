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

const DeleteButton = ({
  isRevision,
  ...props
}: Partial<IMenuItemProps> & { isRevision: boolean }) => {
  return (
    <MenuItem
      {...props}
      shouldDismissPopover={false}
      icon={IconNames.TRASH}
      intent={Intent.DANGER}
      text={isRevision ? "Cancel Revision" : "Delete Collection"}
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
  isRevision: boolean;
}

const StyledMenu = styled(RawMenu)`
  border-radius: 3px;
  padding: 8px;
  & > li:last-child {
    margin-bottom: 0;
  }
`;

const Menu = ({ id = "", isRevision }: Props) => {
  return (
    <StyledMenu>
      <CreateCollection id={id} Button={EditButton} />
      <DeleteCollection id={id} isRevision={isRevision} Button={DeleteButton} />
    </StyledMenu>
  );
};

export default Menu;
