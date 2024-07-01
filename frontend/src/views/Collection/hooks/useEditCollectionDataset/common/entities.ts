import {
  PathParams,
  SubmitOptions,
} from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/hooks/types";

export interface FieldValues {
  [key: string]: string;
}

export type OnEditDatasetFn = (
  pathParams: PathParams,
  fieldValues: FieldValues,
  submitOptions?: SubmitOptions
) => Promise<void>;

export interface EditDatasetAction {
  onEditDataset: OnEditDatasetFn;
}

export interface UseEditCollectionDataset {
  editDatasetAction: EditDatasetAction;
}
