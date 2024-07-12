import { Collection } from "src/common/entities";
import DeleteCollection from "src/components/Collections/components/DeleteCollection";
import CreateCollection from "src/components/CreateCollectionModal";
import { DeleteCollectionFn } from "src/views/Collection/components/CollectionActions";
import { Icon, MenuProps as SDSMenuProps } from "@czi-sds/components";
import {
  StyledMenu,
  StyledMenuItem,
  StyledMenuItemProps,
} from "src/views/Collection/components/CollectionActions/components/MoreDropdown/components/Menu/style";
import {
  DELETE_ICON_PROPS,
  EDIT_ICON_PROPS,
  MENU_PROPS,
} from "src/views/Collection/components/CollectionActions/components/MoreDropdown/components/Menu/constants";
import { Reorder } from "src/views/Collection/hooks/useReorder/common/entities";
import { MENU_ITEM_COLOR } from "src/views/Collection/components/CollectionActions/components/MoreDropdown/components/Menu/types";
import IconSort from "src/views/Collection/components/CollectionActions/components/MoreDropdown/components/Menu/components/IconSort";

interface MenuProps extends Partial<Omit<SDSMenuProps, "onClose">> {
  onClose: () => void;
}

const DeleteButton = ({
  isRevision,
  ...props
}: Partial<StyledMenuItemProps<"TrashCan">> & { isRevision: boolean }) => {
  return (
    <StyledMenuItem
      {...props}
      color={MENU_ITEM_COLOR.ERROR} // Targets custom menu item text color.
      data-testid={
        isRevision ? "dropdown-cancel-revision" : "dropdown-delete-collection"
      }
      icon={<Icon {...DELETE_ICON_PROPS} />}
    >
      {isRevision ? "Cancel Revision" : "Delete Collection"}
    </StyledMenuItem>
  );
};

const EditButton = (props: Partial<StyledMenuItemProps<"Edit">>) => {
  return (
    <StyledMenuItem
      {...props}
      data-testid="dropdown-edit-details"
      icon={<Icon {...EDIT_ICON_PROPS} />}
    >
      Edit Details
    </StyledMenuItem>
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
    <StyledMenu {...MENU_PROPS} {...menuProps} open={Boolean(menuProps.open)}>
      <CreateCollection id={collection.id} Button={EditButton} />
      <StyledMenuItem
        color={MENU_ITEM_COLOR.GRAY} // Targets menu item custom icon color.
        data-testid="dropdown-reorder-datasets"
        disabled={reorder.disabled}
        icon={<IconSort />}
        onClick={() => {
          menuProps.onClose();
          reorder.startReorder();
        }}
      >
        Reorder Datasets
      </StyledMenuItem>
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
