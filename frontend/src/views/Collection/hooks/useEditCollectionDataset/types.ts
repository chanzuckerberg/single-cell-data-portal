import {
  CollectionDatasetFormMethod,
  FieldValues,
  PathParams,
  SubmitOptions,
} from "src/views/Collection/hooks/useCollectionDatasetForm/types";

export interface EditDataset {
  formMethod: CollectionDatasetFormMethod;
  onEditDataset: OnEditDatasetFn;
}

export type OnEditDatasetFn = (
  pathParams: PathParams,
  fieldValues: FieldValues,
  submitOptions?: SubmitOptions
) => Promise<void>;

export interface EditCollectionDatasetAction {
  onEditDataset: OnEditDatasetFn;
}

export interface UseEditCollectionDataset {
  editDatasetAction: EditCollectionDatasetAction;
  formMethod: CollectionDatasetFormMethod;
}
