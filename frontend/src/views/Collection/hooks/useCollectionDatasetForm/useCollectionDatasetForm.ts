import { FormEvent, useCallback, useState } from "react";
import {
  FieldErrors,
  FieldValues,
  PathParams,
  SubmitOptions,
  UseCollectionDatasetForm,
} from "src/views/Collection/hooks/useCollectionDatasetForm/types";
import {
  getAndValidateFieldValues,
  isValid,
} from "src/views/Collection/hooks/useCollectionDatasetForm/utils";

/**
 * Collection dataset form functionality for collection dataset.
 * @returns collection dataset form functionality.
 */
export function useCollectionDatasetForm(): UseCollectionDatasetForm {
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
    []
  );

  return { clearErrors, errors, handleSubmit };
}
