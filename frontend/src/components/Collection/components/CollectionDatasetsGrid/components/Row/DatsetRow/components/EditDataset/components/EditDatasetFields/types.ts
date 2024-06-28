import { FieldValues } from "src/views/Collection/hooks/useEditCollectionDataset/common/entities";
import { FieldErrors } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/EditDataset/hooks/types";

export interface Props {
  errors: FieldErrors;
  fieldValues: FieldValues;
}
