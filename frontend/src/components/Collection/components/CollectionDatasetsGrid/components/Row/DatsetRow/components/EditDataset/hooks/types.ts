import { FormEvent } from "react";
import { FieldValues } from "src/views/Collection/hooks/useEditCollectionDataset/common/entities";
import { Collection, Dataset } from "src/common/entities";

export type FieldErrors = Record<keyof FieldValues, string>;

export type HandleSubmit = (
  onSubmit: OnSubmit,
  pathParams: PathParams,
  defaultValues: FieldValues,
  submitOptions?: SubmitOptions
) => (event: FormEvent) => Promise<void>;

export type OnSubmit = (
  pathParams: PathParams,
  fieldValues: FieldValues,
  submitOptions?: SubmitOptions
) => Promise<void>;

export interface PathParams {
  collectionId: Collection["id"];
  datasetId: Dataset["id"];
}

export interface SubmitOptions {
  onError?: () => void;
  onSuccess?: () => void;
}
