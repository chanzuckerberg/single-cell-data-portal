import React, { useEffect, useMemo } from "react";
import {
  FIELD_NAMES,
  INPUT_TEXT_PROPS,
} from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/components/EditDatasetFields/constants";
import { Props } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/components/EditDatasetFields/types";
import {
  StyledInputText,
  StyledLabel,
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
      {...INPUT_TEXT_PROPS.TITLE}
      defaultValue={defaultValues[FIELD_NAMES.TITLE]}
      error={Boolean(errors[FIELD_NAMES.TITLE])}
      helperText={errors[FIELD_NAMES.TITLE]}
      label={<StyledLabel>Name</StyledLabel>}
      onChange={(e) => {
        register(FIELD_NAMES.TITLE).onChange(e);
        if (errors[FIELD_NAMES.TITLE]) clearErrors(FIELD_NAMES.TITLE);
      }}
      onKeyDown={(e) => e.stopPropagation()} // Prevents the input field from losing focus when a key matching the first letter of a menu item is pressed.
    />
  );
}
