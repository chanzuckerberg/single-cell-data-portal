import {
  ChangeEvent,
  FormEvent,
  useCallback,
  useEffect,
  useState,
} from "react";
import {
  DefaultValues,
  DirtyFields,
  FieldErrors,
  FieldName,
  FieldValue,
  FieldValues,
  PathParams,
  SubmitOptions,
  UseCollectionDatasetForm,
} from "src/views/Collection/hooks/useCollectionDatasetForm/types";
import {
  getAndValidateFieldValues,
  isValid,
  updateDirtyFields,
  updateFieldErrors,
} from "src/views/Collection/hooks/useCollectionDatasetForm/utils";
import {
  DEFAULT_DIRTY_FIELDS,
  DEFAULT_FIELD_ERRORS,
  DEFAULT_FIELD_VALUES,
} from "src/views/Collection/hooks/useCollectionDatasetForm/constants";

/**
 * Collection dataset form functionality for collection dataset.
 * @returns collection dataset form functionality.
 */
export function useCollectionDatasetForm(): UseCollectionDatasetForm {
  const [defaultValues, setDefaultValues] =
    useState<DefaultValues>(DEFAULT_FIELD_VALUES);
  const [dirtyFields, setDirtyFields] =
    useState<DirtyFields>(DEFAULT_DIRTY_FIELDS);
  const [errors, setErrors] = useState<FieldErrors>(DEFAULT_FIELD_ERRORS);
  const [isDirty, setIsDirty] = useState<boolean>(false);

  const clearErrors = useCallback((name?: FieldName) => {
    setErrors((prevErrors) => updateFieldErrors(prevErrors, name));
  }, []);

  const handleSubmit = useCallback(
    (
      onSubmit: (
        pathParams: PathParams,
        fieldValues: FieldValues,
        submitOptions?: SubmitOptions
      ) => Promise<void>,
      pathParams: PathParams,
      submitOptions?: SubmitOptions
    ) => {
      return async (event: FormEvent) => {
        event.preventDefault();
        const [fieldValues, fieldErrors] = getAndValidateFieldValues(
          event.target as HTMLFormElement,
          defaultValues
        );
        if (isValid(fieldErrors)) {
          await onSubmit(pathParams, fieldValues, submitOptions);
        } else {
          setErrors(fieldErrors);
        }
      };
    },
    [defaultValues]
  );

  const register = useCallback(
    (name: FieldName) => {
      return {
        onChange: (
          event: ChangeEvent<HTMLInputElement | HTMLTextAreaElement>
        ) => {
          const value = event.target.value as FieldValue;
          setDirtyFields((prevDirtyFields) =>
            updateDirtyFields(prevDirtyFields, { [name]: value }, defaultValues)
          );
        },
      };
    },
    [defaultValues]
  );

  useEffect(() => {
    setIsDirty(Object.keys(dirtyFields).length > 0);
  }, [dirtyFields]);

  return {
    clearErrors,
    formState: { dirtyFields, errors, isDirty },
    handleSubmit,
    register,
    setDefaultValues,
  };
}
