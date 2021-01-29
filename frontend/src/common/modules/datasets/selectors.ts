import { Dataset } from "src/common/entities";

export function hasAssets(dataset?: Dataset) {
  if (!dataset) return false;

  return Boolean(dataset?.dataset_assets?.length);
}
