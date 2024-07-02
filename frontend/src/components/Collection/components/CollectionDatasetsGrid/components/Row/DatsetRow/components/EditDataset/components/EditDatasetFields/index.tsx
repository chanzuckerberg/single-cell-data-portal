import React, { useEffect, useMemo } from "react";
import {
  TITLE_FIELD_NAME,
  TITLE_INPUT_TEXT_PROPS,
} from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/components/EditDatasetFields/constants";
import { Props } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/components/EditDatasetFields/types";
import {
  Label,
  StyledInputText,
} from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/components/EditDatasetFields/style";
import { mapDatasetToFormDefaultValues } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/components/EditDatasetFields/utils";

export default function EditDatasetFields({
  dataset,
  formMethod,
}: Props): JSX.Element {
  const {
    clearErrors,
    formState: { errors },
    register,
    setDefaultValues,
  } = formMethod;
  const defaultValues = useMemo(
    () => mapDatasetToFormDefaultValues(dataset),
    [dataset]
  );

  useEffect(() => {
    setDefaultValues(defaultValues);
  }, [defaultValues, setDefaultValues]);

  return (
    <StyledInputText
      {...TITLE_INPUT_TEXT_PROPS}
      defaultValue={defaultValues.title}
      error={Boolean(errors.title)}
      helperText={errors.title}
      label={<Label>Name</Label>}
      onChange={(e) => {
        register(TITLE_FIELD_NAME).onChange(e);
        if (errors.title) clearErrors(TITLE_FIELD_NAME);
      }}
      onKeyDown={(e) => e.stopPropagation()} // Prevents the input field from losing focus when a key matching the first letter of a menu item is pressed.
    />
  );
}
