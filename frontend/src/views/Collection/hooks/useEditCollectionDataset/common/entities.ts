import { Dataset } from "src/common/entities";

export interface FieldValues {
  [key: string]: string;
}

export type OnEditDatasetFn = (
  dataset: Dataset,
  fieldValues: FieldValues
) => Promise<void>;

export interface EditDatasetAction {
  onEditDataset: OnEditDatasetFn;
}

export interface UseEditCollectionDataset {
  editDatasetAction: EditDatasetAction;
}
