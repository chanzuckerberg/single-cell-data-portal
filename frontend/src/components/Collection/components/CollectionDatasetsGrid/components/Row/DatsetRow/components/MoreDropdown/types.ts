import { Collection, Dataset } from "src/common/entities";
import { MenuItemProps } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/MoreDropdown/components/Menu/types";

export interface Props {
  collectionId: Collection["id"];
  dataset: Dataset;
  disabled: boolean;
  menuItemProps: MenuItemProps;
}
