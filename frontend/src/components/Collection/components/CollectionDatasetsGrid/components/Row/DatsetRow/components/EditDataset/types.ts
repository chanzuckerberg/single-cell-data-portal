import { ElementType } from "react";
import { Collection, Dataset } from "src/common/entities";
import { UseMenu } from "src/views/Collection/hooks/useMenu";
import { EditDataset } from "src/views/Collection/hooks/useEditCollectionDataset/types";

export interface Props {
  Button: ElementType;
  collectionId: Collection["id"];
  dataset: Dataset;
  editDataset: EditDataset;
  menuProps: Omit<UseMenu, "onOpen">;
}
