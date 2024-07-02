import { FieldValues } from "src/views/Collection/hooks/useCollectionDatasetForm/types";
import { Dataset } from "src/common/entities";

/**
 * Maps dataset to form default values.
 * @param dataset - Dataset.
 * @returns dataset mapped to form default values.
 */
export function mapDatasetToFormDefaultValues(dataset: Dataset): FieldValues {
  return {
    title: dataset.name,
  };
}
