import { FormEvent, useCallback, useState } from "react";
import {
  FieldErrors,
  HandleSubmit,
} from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/hooks/types";
import { Dataset } from "src/common/entities";
import { FieldValues } from "src/views/Collection/hooks/useEditCollectionDataset/common/entities";

export interface UseForm {
  clearErrors: () => void;
  errors: FieldErrors;
  handleSubmit: HandleSubmit;
}

/**
 * Edit dataset form functionality for collection datasets.
 * @returns edit dataset form functionality.
 */
export function useForm(): UseForm {
  const [errors, setErrors] = useState<FieldErrors>({} as FieldErrors);

  const clearErrors = useCallback(() => {
    setErrors({} as FieldErrors);
  }, []);

  const handleSubmit = useCallback(
    (
      onSubmit: (dataset: Dataset, fieldValues: FieldValues) => Promise<void>,
      dataset: Dataset,
      defaultValues: FieldValues
    ) => {
      return async (event: FormEvent) => {
        event.preventDefault();
        const fieldValues = getFieldValues(
          event.target as HTMLFormElement,
          defaultValues
        );
        const errors = validateForm(fieldValues, defaultValues);
        if (isValid(errors)) {
          onSubmit(dataset, fieldValues);
        } else {
          setErrors(errors);
        }
      };
    },
    []
  );

  return { clearErrors, errors, handleSubmit };
}

/**
 * Returns field values from form data.
 * @param formEl - Form element.
 * @param defaultValues - Default field values.
 * @returns field values.
 */
function getFieldValues(
  formEl: HTMLFormElement,
  defaultValues: FieldValues
): FieldValues {
  const formData: FormData = new FormData(formEl);
  const fieldValues = {} as FieldValues;
  for (const key of Object.keys(defaultValues)) {
    if (formData.has(key)) {
      Object.assign(fieldValues, { [key]: formData.get(key) });
    }
  }
  return fieldValues;
}

/**
 * Returns true if the form data is valid.
 * @param errors - Field errors.
 * @returns true if the form data is valid.
 */
function isValid(errors: FieldErrors): boolean {
  return Object.keys(errors).length === 0;
}

/**
 * Validates form data.
 * @param fieldValues - Field values.
 * @param defaultValues - Default field values.
 * @returns error.
 */
function validateForm(
  fieldValues: FieldValues,
  defaultValues: FieldValues
): FieldErrors {
  const errors = {} as FieldErrors;
  for (const [key, value] of Object.entries(fieldValues)) {
    if (!value || value === defaultValues[key]) {
      Object.assign(errors, { [key]: `${key} is required` });
    }
  }
  return errors;
}
