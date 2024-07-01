import { FormEvent, useCallback, useState } from "react";
import {
  FieldErrors,
  HandleSubmit,
  PathParams,
  SubmitOptions,
} from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/hooks/types";
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
      onSubmit: (
        pathParams: PathParams,
        fieldValues: FieldValues,
        submitOptions?: SubmitOptions
      ) => Promise<void>,
      pathParams: PathParams,
      defaultValues: FieldValues,
      submitOptions?: SubmitOptions
    ) => {
      return async (event: FormEvent) => {
        event.preventDefault();
        const fieldValues = getFieldValues(
          event.target as HTMLFormElement,
          defaultValues
        );
        const errors = validateForm(fieldValues, defaultValues);
        if (isValid(errors)) {
          await onSubmit(pathParams, fieldValues, submitOptions);
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
      const value = formData.get(key);
      Object.assign(fieldValues, { [key]: formatFormDataValue(value) });
    }
  }
  return fieldValues;
}

/**
 * Returns form data entry value, formatted, if necessary.
 * @param value - Form data entry value.
 * @returns formatted data entry value.
 */
function formatFormDataValue(
  value: FormDataEntryValue | null
): FormDataEntryValue | null {
  if (typeof value === "string") {
    return value.trim();
  }
  return value;
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
