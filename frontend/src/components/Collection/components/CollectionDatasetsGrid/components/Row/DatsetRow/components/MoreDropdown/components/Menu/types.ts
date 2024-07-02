import { Props as ChooserProps } from "src/components/DropboxChooser";
import { Collection, Dataset } from "src/common/entities";
import { UseMenu } from "src/views/Collection/hooks/useMenu";
import { EditDataset } from "src/views/Collection/hooks/useEditCollectionDataset/types";

export interface MenuItemProps {
  editDataset: EditDataset;
  isFailed: boolean;
  isLoading: boolean;
  onUploadFile: ChooserProps["onUploadFile"];
  revisionsEnabled: boolean;
}

export interface Props {
  collectionId: Collection["id"];
  dataset: Dataset;
  menuItemProps: MenuItemProps;
  menuProps: Omit<UseMenu, "onOpen">;
}
