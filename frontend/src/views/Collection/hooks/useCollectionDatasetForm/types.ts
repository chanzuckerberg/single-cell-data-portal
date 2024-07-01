import { FormEvent } from "react";

export type FieldErrors = Record<keyof FieldValues, string>;

export interface FieldValues {
  [key: string]: string;
}

export type CollectionDatasetFormMethod = UseCollectionDatasetForm;

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

export type PathParams = Record<string, string>;

export interface SubmitOptions {
  onError?: () => void;
  onSuccess?: () => void;
}

export interface UseCollectionDatasetForm {
  clearErrors: () => void;
  errors: FieldErrors;
  handleSubmit: HandleSubmit;
}
