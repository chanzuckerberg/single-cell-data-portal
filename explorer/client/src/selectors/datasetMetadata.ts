/**
 * Selectors for accessing and deriving values from datasetMetadata slice of store.
 */

// App dependencies
import { RootState } from "../reducers";
import { DatasetMetadataState } from "../reducers/datasetMetadata";

export const selectDatasetMetadata = (state: RootState): DatasetMetadataState =>
  state.datasetMetadata;

/**
 * Determine if seamless features are enabled.
 * @param state - root state of store.
 * @return True if both dataset metadata and Portal URL are defined.
 */
export const selectIsSeamlessEnabled = (state: RootState): boolean => {
  const datasetMetadata = selectDatasetMetadata(state);
  return Boolean(datasetMetadata.datasetMetadata && datasetMetadata.portalUrl);
};
