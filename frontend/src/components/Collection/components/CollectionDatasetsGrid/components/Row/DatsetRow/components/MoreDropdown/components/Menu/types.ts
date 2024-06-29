import { Props as ChooserProps } from "src/components/DropboxChooser";
import { Collection, Dataset } from "src/common/entities";
import { UseMoreMenu } from "src/views/Collection/hooks/useMoreMenu";

export interface MenuItemProps {
  isLoading: boolean;
  onUploadFile: ChooserProps["onUploadFile"];
  revisionsEnabled: boolean;
}

export interface Props {
  collectionId: Collection["id"];
  dataset: Dataset;
  menuItemProps: MenuItemProps;
  menuProps: Omit<UseMoreMenu, "onOpen">;
}
