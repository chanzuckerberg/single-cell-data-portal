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
  Menu as StyledMenu,
  MenuItem,
  ReorderMenuItem,
} from "src/views/Collection/components/CollectionActions/components/MoreDropdown/components/Menu/style";
import { DEFAULT_MENU_PROPS } from "src/views/Collection/components/CollectionActions/components/MoreDropdown/components/Menu/constants";
import IconSort from "src/views/Collection/components/CollectionActions/components/MoreDropdown/components/Menu/components/IconSort";
import { Reorder } from "src/views/Collection/hooks/useReorder/common/entities";

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
    <MenuItem
      {...props}
      data-testid="dropdown-edit-details"
      sdsIcon="edit"
      sdsIconProps={{ color: "gray", shade: 400 }}
    >
      Edit Details
    </MenuItem>
  );
};

interface Props {
  collection: Collection;
  handleDeleteCollection: DeleteCollectionFn;
  isDeleting: boolean;
  isRevision: boolean;
  menuProps: MenuProps;
  reorder: Reorder;
}

const Menu = ({
  collection,
  handleDeleteCollection,
  isDeleting,
  isRevision,
  menuProps,
  reorder,
}: Props) => {
  return (
    <StyledMenu
      {...DEFAULT_MENU_PROPS}
      {...menuProps}
      open={Boolean(menuProps.open)}
    >
      <CreateCollection id={collection.id} Button={EditButton} />
      {reorder.isReorderUX && (
        <ReorderMenuItem
          disabled={reorder.disabled}
          onClick={() => {
            menuProps.onClose();
            reorder.startReorder();
          }}
        >
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
