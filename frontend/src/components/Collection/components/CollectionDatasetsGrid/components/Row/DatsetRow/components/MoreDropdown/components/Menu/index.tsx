import DropboxChooser from "src/components/DropboxChooser";
import { Props } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/MoreDropdown/components/Menu/types";
import DeleteDataset from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/DeleteDataset";
import {
  StyledMenu,
  StyledMenuItem,
} from "src/views/Collection/components/CollectionActions/components/MoreDropdown/components/Menu/style";
import {
  DELETE_ICON_PROPS,
  EDIT_ICON_PROPS,
  MENU_PROPS,
} from "src/views/Collection/components/CollectionActions/components/MoreDropdown/components/Menu/constants";
import { MENU_ITEM_COLOR } from "src/views/Collection/components/CollectionActions/components/MoreDropdown/components/Menu/types";
import { Icon } from "@czi-sds/components";
import EditDataset from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset";
import { isDeleteDatasetAvailable } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/MoreDropdown/components/Menu/utils";
import { DROPDOWN_EDIT_DATASET } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/MoreDropdown/components/Menu/constants";

export default function Menu({
  collectionId,
  dataset,
  menuItemProps,
  menuProps,
}: Props): JSX.Element {
  const { editDataset, isFailed, isLoading, onUploadFile, revisionsEnabled } =
    menuItemProps;
  const canDelete = isDeleteDatasetAvailable(dataset, revisionsEnabled);
  return (
    <StyledMenu {...MENU_PROPS} {...menuProps}>
      {/* Edit */}
      <EditDataset
        Button={(buttonProps) => (
          <StyledMenuItem
            data-testid={DROPDOWN_EDIT_DATASET}
            disabled={isFailed || isLoading}
            icon={<Icon {...EDIT_ICON_PROPS} />}
            {...buttonProps}
          >
            Rename Dataset
          </StyledMenuItem>
        )}
        collectionId={collectionId}
        dataset={dataset}
        editDataset={editDataset}
        menuProps={menuProps}
      />
      {/* Upload */}
      {revisionsEnabled && (
        <DropboxChooser onUploadFile={onUploadFile}>
          <StyledMenuItem
            disabled={isLoading}
            icon={<Icon {...EDIT_ICON_PROPS} />}
          >
            Update Dataset
          </StyledMenuItem>
        </DropboxChooser>
      )}
      {/* Delete */}
      {canDelete && (
        <DeleteDataset
          Button={(buttonProps) => (
            <StyledMenuItem
              color={MENU_ITEM_COLOR.ERROR} // Targets custom menu item text color.
              icon={<Icon {...DELETE_ICON_PROPS} />}
              {...buttonProps}
            >
              Delete Dataset
            </StyledMenuItem>
          )}
          collectionId={collectionId}
          datasetId={dataset.id}
        />
      )}
    </StyledMenu>
  );
}
