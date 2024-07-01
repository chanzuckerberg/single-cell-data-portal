import { ElementType } from "react";
import { Collection, Dataset } from "src/common/entities";
import { UseMoreMenu } from "src/views/Collection/hooks/useMoreMenu";
import { EditDataset } from "src/views/Collection/hooks/useEditCollectionDataset/types";

export interface Props {
  Button: ElementType;
  collectionId: Collection["id"];
  dataset: Dataset;
  editDataset: EditDataset;
  menuProps: Omit<UseMoreMenu, "onOpen">;
}
