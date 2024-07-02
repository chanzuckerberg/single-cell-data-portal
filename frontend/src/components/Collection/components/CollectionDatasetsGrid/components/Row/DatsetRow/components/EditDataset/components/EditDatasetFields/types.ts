import { CollectionDatasetFormMethod } from "src/views/Collection/hooks/useCollectionDatasetForm/types";
import { Dataset } from "src/common/entities";

export interface Props {
  dataset: Dataset;
  formMethod: CollectionDatasetFormMethod;
}
