import { FormEvent } from "react";
import { FieldValues } from "src/views/Collection/hooks/useEditCollectionDataset/common/entities";
import { Dataset } from "src/common/entities";

export type FieldErrors = Record<keyof FieldValues, string>;

export type HandleSubmit = (
  onSubmit: OnSubmit,
  dataset: Dataset,
  defaultValues: FieldValues
) => (event: FormEvent) => Promise<void>;

export type OnSubmit = (
  dataset: Dataset,
  fieldValues: FieldValues
) => Promise<void>;
