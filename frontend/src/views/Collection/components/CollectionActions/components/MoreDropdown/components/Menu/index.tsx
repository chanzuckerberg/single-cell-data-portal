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
import { DeleteCollectionFn } from "src/views/Collection/components/CollectionActions";

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
  collection: Collection;
  handleDeleteCollection: DeleteCollectionFn;
  isDeleting: boolean;
  isRevision: boolean;
}

const StyledMenu = styled(RawMenu)`
  border-radius: 3px;
  padding: 8px;

  & > li:last-child {
    margin-bottom: 0;
  }
`;

const Menu = ({
  collection,
  handleDeleteCollection,
  isDeleting,
  isRevision,
}: Props) => {
  return (
    <StyledMenu>
      <CreateCollection id={collection.id} Button={EditButton} />
      <DeleteCollection
        Button={DeleteButton}
        handleDeleteCollection={handleDeleteCollection}
        isDeleting={isDeleting}
        isRevision={isRevision}
      />
    </StyledMenu>
  );
};

export default Menu;
