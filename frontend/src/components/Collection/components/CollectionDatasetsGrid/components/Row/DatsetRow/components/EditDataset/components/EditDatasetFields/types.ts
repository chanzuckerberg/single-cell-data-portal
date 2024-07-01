import {
  CollectionDatasetFormMethod,
  FieldValues,
} from "src/views/Collection/hooks/useCollectionDatasetForm/types";

export interface Props {
  defaultValues: FieldValues;
  formMethod: CollectionDatasetFormMethod;
}
