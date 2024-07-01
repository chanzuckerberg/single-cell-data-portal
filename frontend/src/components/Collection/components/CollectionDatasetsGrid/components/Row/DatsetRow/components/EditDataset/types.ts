import { ElementType } from "react";
import { Collection, Dataset } from "src/common/entities";
import { UseMoreMenu } from "src/views/Collection/hooks/useMoreMenu";

export interface Props {
  Button: ElementType;
  collectionId: Collection["id"];
  dataset: Dataset;
  menuProps: Omit<UseMoreMenu, "onOpen">;
}
