import { ElementType } from "react";
import { Collection, Dataset } from "src/common/entities";

export interface Props {
  Button: ElementType;
  collectionId: Collection["id"];
  dataset: Dataset;
}
