import { ChangeEvent, Dispatch, FormEvent, SetStateAction } from "react";

export type CollectionDatasetFormMethod = UseCollectionDatasetForm;

export type DefaultValues = Partial<FieldValues>;

export type DirtyFields = Partial<Record<FieldName, boolean>>;

export type FieldErrors = Partial<Record<FieldName, string>>;

export type FieldName = keyof FieldValues;

export type FieldValue = string; // Currently, only string values are supported.

export type FieldValues = Record<string, FieldValue>;

export interface FormState {
  dirtyFields: DirtyFields;
  errors: FieldErrors;
  isDirty: boolean;
}

export type HandleSubmit = (
  onSubmit: OnSubmit,
  pathParams: PathParams,
  submitOptions?: SubmitOptions
) => (event: FormEvent) => Promise<void>;

export type OnSubmit = (
  pathParams: PathParams,
  fieldValues: FieldValues,
  submitOptions?: SubmitOptions
) => Promise<void>;

export type PathParams = Record<string, string>;

export type RegisterFn = (name: FieldName) => {
  onChange: (
    event: ChangeEvent<HTMLInputElement | HTMLTextAreaElement>
  ) => void;
};

export interface SubmitOptions {
  onError?: () => void;
  onSuccess?: () => void;
}

export interface UseCollectionDatasetForm {
  clearErrors: (name?: FieldName) => void;
  formState: FormState;
  handleSubmit: HandleSubmit;
  register: RegisterFn;
  reset: () => void;
  setDefaultValues: Dispatch<SetStateAction<DefaultValues>>;
}
