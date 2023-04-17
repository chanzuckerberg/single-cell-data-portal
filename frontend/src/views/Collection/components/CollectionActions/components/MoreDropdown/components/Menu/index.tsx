import {
  IMenuItemProps,
  Intent,
  Menu as RawMenu,
  MenuItem,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import styled from "@emotion/styled";
import { Collection } from "src/common/entities";
import DeleteCollection from "src/components/Collections/components/DeleteCollection";
import CreateCollection from "src/components/CreateCollectionModal";

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
      data-testid={
        isRevision ? "dropdown-cancel-revision" : "dropdown-delete-collection"
      }
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
      data-testid="dropdown-edit-details"
    />
  );
};

interface Props {
  id?: string;
  isRevision: boolean;
  visibility: Collection["visibility"];
}

const StyledMenu = styled(RawMenu)`
  border-radius: 3px;
  padding: 8px;
  & > li:last-child {
    margin-bottom: 0;
  }
`;

const Menu = ({ id = "", isRevision, visibility }: Props) => {
  return (
    <StyledMenu>
      <CreateCollection id={id} Button={EditButton} />
      <DeleteCollection
        id={id}
        isRevision={isRevision}
        Button={DeleteButton}
        visibility={visibility}
      />
    </StyledMenu>
  );
};

export default Menu;
