import { Collection } from "src/common/entities";
import DeleteCollection from "src/components/Collections/components/DeleteCollection";
import CreateCollection from "src/components/CreateCollectionModal";
import { DeleteCollectionFn } from "src/views/Collection/components/CollectionActions";
import {
  MenuItemProps as SDSMenuItemProps,
  MenuProps as SDSMenuProps,
} from "@czi-sds/components";
import {
  DeleteMenuItem,
  EditMenuItem,
  Menu as StyledMenu,
  ReorderMenuItem,
} from "src/views/Collection/components/CollectionActions/components/MoreDropdown/components/Menu/style";
import { DEFAULT_MENU_PROPS } from "src/views/Collection/components/CollectionActions/components/MoreDropdown/components/Menu/constants";
import IconSort from "src/views/Collection/components/CollectionActions/components/MoreDropdown/components/Menu/components/IconSort";
import { OnReorderFn, REORDER_MODE } from "src/common/hooks/useReorderMode";

interface MenuProps extends Partial<Omit<SDSMenuProps, "onClose">> {
  onClose: () => void;
}

const DeleteButton = ({
  isRevision,
  ...props
}: Partial<SDSMenuItemProps<"trashCan">> & { isRevision: boolean }) => {
  return (
    <DeleteMenuItem
      {...props}
      data-testid={
        isRevision ? "dropdown-cancel-revision" : "dropdown-delete-collection"
      }
      sdsIcon="trashCan"
      sdsIconProps={{ color: "error" }}
    >
      {isRevision ? "Cancel Revision" : "Delete Collection"}
    </DeleteMenuItem>
  );
};

const EditButton = (props: Partial<SDSMenuItemProps<"edit">>) => {
  return (
    <EditMenuItem
      {...props}
      data-testid="dropdown-edit-details"
      sdsIcon="edit"
      sdsIconProps={{ color: "gray" }}
    >
      Edit Details
    </EditMenuItem>
  );
};

interface Props {
  collection: Collection;
  handleDeleteCollection: DeleteCollectionFn;
  isDeleting: boolean;
  isReorderUX: boolean;
  isRevision: boolean;
  menuProps: MenuProps;
  onReorder: OnReorderFn;
}

const Menu = ({
  collection,
  handleDeleteCollection,
  isDeleting,
  isReorderUX,
  isRevision,
  menuProps,
  onReorder,
}: Props) => {
  // Facilitates dataset reordering functionality.
  const handleReorder = () => {
    menuProps.onClose();
    onReorder(REORDER_MODE.ACTIVE);
  };
  return (
    <StyledMenu
      {...DEFAULT_MENU_PROPS}
      {...menuProps}
      open={Boolean(menuProps.open)}
    >
      <CreateCollection id={collection.id} Button={EditButton} />
      {isReorderUX && (
        <ReorderMenuItem onClick={handleReorder}>
          <IconSort />
          Reorder Datasets
        </ReorderMenuItem>
      )}
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
