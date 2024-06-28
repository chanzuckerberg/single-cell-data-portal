import { Props as ChooserProps } from "src/components/DropboxChooser";
import { Dataset } from "src/common/entities";
import { UseMoreMenu } from "src/views/Collection/hooks/useMoreMenu";

export interface MenuItemProps {
  isLoading: boolean;
  onUploadFile: ChooserProps["onUploadFile"];
  revisionsEnabled: boolean;
}

export interface Props {
  dataset: Dataset;
  menuItemProps: MenuItemProps;
  menuProps: Omit<UseMoreMenu, "onOpen">;
}
