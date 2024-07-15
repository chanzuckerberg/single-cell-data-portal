import {
  DefaultValues,
  DirtyFields,
  FieldErrors,
  FieldValues,
} from "src/views/Collection/hooks/useCollectionDatasetForm/types";
import { DEFAULT_FIELD_ERRORS } from "src/views/Collection/hooks/useCollectionDatasetForm/constants";

/**
 * Clears given field error or all field errors.
 * @param fieldErrors - Field errors.
 * @param name - Field name.
 * @returns updated field errors.
 */
export function clearFieldErrors(
  fieldErrors: FieldErrors,
  name?: string
): FieldErrors {
  const updatedFieldErrors = { ...fieldErrors };
  if (!name) return DEFAULT_FIELD_ERRORS;
  if (name in updatedFieldErrors) {
    delete updatedFieldErrors[name];
  }
  return updatedFieldErrors;
}

/**
 * Returns field values and errors from form data.
 * @param formEl - Form element.
 * @returns field values and errors tuple.
 */
export function getAndValidateFieldValues(
  formEl: HTMLFormElement
): [FieldValues, FieldErrors] {
  const fieldValues = getFieldValues(formEl);
  const fieldErrors = validateForm(formEl, fieldValues);
  return [fieldValues, fieldErrors];
}

/**
 * Returns field label value.
 * @param formEl - Form element.
 * @param key - Field value key.
 * @returns field label value.
 */
function getFieldLabelValue(formEl: HTMLFormElement, key: string): string {
  const labelEl = formEl.querySelector(`label[for="${key}"]`);
  return labelEl?.textContent || key;
}

/**
 * Returns field values from form data.
 * @param formEl - Form element.
 * @returns field values.
 */
export function getFieldValues(formEl: HTMLFormElement): FieldValues {
  const formData: FormData = new FormData(formEl);
  const fieldValues = {} as FieldValues;
  formData.forEach((value: FormDataEntryValue, key: string) => {
    // Currently, only string values are supported.
    if (typeof value === "string") {
      Object.assign(fieldValues, { [key]: formatFormDataValue(value) });
    }
  });
  return fieldValues;
}

/**
 * Returns form data entry value, formatted, if necessary.
 * @param value - Form data entry value.
 * @returns formatted data entry value.
 */
export function formatFormDataValue(value: string): string {
  return value.trim();
}

/**
 * Returns true if the form data is valid.
 * @param errors - Field errors.
 * @returns true if the form data is valid.
 */
export function isValid(errors: FieldErrors): boolean {
  return Object.keys(errors).length === 0;
}

/**
 * Updates dirty fields.
 * If the field value is different from the default value, the field is marked as dirty.
 * @param dirtyFields - Dirty fields.
 * @param fieldValues - Field values.
 * @param defaultValues - Default field values.
 * @returns updated dirty fields.
 */
export function updateDirtyFields(
  dirtyFields: DirtyFields,
  fieldValues: FieldValues,
  defaultValues: DefaultValues
): DirtyFields {
  const updatedDirtyFields = { ...dirtyFields };
  for (const [key, value] of Object.entries(fieldValues)) {
    if (value !== defaultValues[key]) {
      Object.assign(updatedDirtyFields, { [key]: true });
    } else {
      delete updatedDirtyFields[key];
    }
  }
  return updatedDirtyFields;
}

/**
 * Validates form data.
 * At present, form validation only checks if the field value is empty.
 * @param formEl - Form element.
 * @param fieldValues - Field values.
 * @returns error.
 */
export function validateForm(
  formEl: HTMLFormElement,
  fieldValues: FieldValues
): FieldErrors {
  const fieldErrors = {} as FieldErrors;
  for (const [key, value] of Object.entries(fieldValues)) {
    if (!value) {
      Object.assign(fieldErrors, {
        [key]: `${getFieldLabelValue(formEl, key)} is required`,
      });
    }
  }
  return fieldErrors;
}
