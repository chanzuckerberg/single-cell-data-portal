import React from "react";
import { TITLE_INPUT_TEXT_PROPS } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/components/EditDatasetFields/constants";
import { Props } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/components/EditDatasetFields/types";
import {
  Label,
  StyledInputText,
} from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/components/EditDatasetFields/style";

export default function EditDatasetFields({
  clearErrors,
  errors,
  fieldValues,
}: Props): JSX.Element {
  return (
    <StyledInputText
      {...TITLE_INPUT_TEXT_PROPS}
      defaultValue={fieldValues.title}
      error={Boolean(errors.title)}
      helperText={errors.title}
      label={<Label>Name</Label>}
      onChange={() => {
        if (errors.title) clearErrors();
      }}
      onKeyDown={(e) => e.stopPropagation()} // Prevents the input field from losing focus when a key matching the first letter of a menu item is pressed.
    />
  );
}
