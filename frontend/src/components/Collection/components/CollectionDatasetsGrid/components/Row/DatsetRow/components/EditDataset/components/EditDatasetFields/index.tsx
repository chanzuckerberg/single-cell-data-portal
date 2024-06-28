import React from "react";
import { TITLE_INPUT_TEXT_PROPS } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/components/EditDatasetFields/constants";
import { Props } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/components/EditDatasetFields/types";
import {
  Label,
  StyledInputText,
} from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/components/EditDatasetFields/style";

export default function EditDatasetFields({
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
    />
  );
}
