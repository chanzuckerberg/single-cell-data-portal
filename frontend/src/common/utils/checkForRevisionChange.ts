import isEmpty from "lodash/isEmpty";
import isEqual from "lodash/isEqual";
import xorWith from "lodash/xorWith";
import { Collection, Dataset } from "../entities";

const IGNORED_COLLECTION_FIELDS = [
  "id",
  "visibility",
  "created_at",
  "updated_at",
  "revisioning_in",
  "revision_of",
  "revision_diff",
  "datasets",
  "genesets",
  "links",
] as Array<keyof Collection>;
const IGNORED_DATASET_FIELDS = [
  "created_at",
  "updated_at",
  "collection_visibility",
  "original_uuid",
  "id",
  "collection_id",
  "processing_status",
  "dataset_assets",
  "dataset_deployments",
] as Array<keyof Dataset>;

function checkListForChanges(
  revisedList: Array<unknown>,
  publishedList: Array<unknown>
): boolean {
  return !isEmpty(xorWith(revisedList, publishedList, isEqual));
}

function checkDatasetKeyForDifference(
  datasetKey: keyof Dataset,
  revisedDataset: Dataset,
  publishedDataset: Dataset
) {
  if (publishedDataset[datasetKey] instanceof Array) {
    // Dataset entry is an array
    if (
      checkListForChanges(
        publishedDataset[datasetKey] as Array<unknown>,
        revisedDataset[datasetKey] as Array<unknown>
      )
    ) {
      console.log(datasetKey);
      return true;
    }
  } else if (
    // Dataset entry is an object
    publishedDataset[datasetKey] instanceof Object
  ) {
    if (!isEqual(publishedDataset[datasetKey], revisedDataset[datasetKey])) {
      console.log(datasetKey);
      return true;
    }
  } else if (
    // scalar value
    publishedDataset[datasetKey] !== revisedDataset[datasetKey]
  ) {
    return true;
  }
  return false;
}

function checkDatasetsForChanges(
  revisedDatasets: Map<Dataset["id"], Dataset>,
  publishedDatasets: Map<Dataset["id"], Dataset>
): boolean {
  if (publishedDatasets.size !== revisedDatasets.size) return true;
  // Check dataset fields for differences
  return Array.from(publishedDatasets.values()).some((publishedDataset) => {
    const revisedDataset = revisedDatasets.get(publishedDataset.id) as Dataset;
    let datasetKey = "" as keyof Dataset;
    for (datasetKey in publishedDataset) {
      if (IGNORED_DATASET_FIELDS.includes(datasetKey)) {
        continue; // ignore these fields (I hate using continue, but sonarjs considers this best)
      }
      if (
        checkDatasetKeyForDifference(
          datasetKey,
          revisedDataset,
          publishedDataset
        )
      )
        return true;
    }

    return false;
  });
}

export default function checkForRevisionChange(
  revision: Collection,
  publishedCollection: Collection
): boolean {
  // Check collection fields for differences
  let collectionKey = "" as keyof Collection;
  for (collectionKey in publishedCollection) {
    if (
      !IGNORED_COLLECTION_FIELDS.includes(collectionKey) &&
      publishedCollection[collectionKey] !== revision[collectionKey]
    ) {
      console.log("found diff for");
      console.log(collectionKey);
      console.log(publishedCollection[collectionKey]);
      console.log(revision[collectionKey]);
      return true;
    }
  }
  if (publishedCollection.links.length !== revision.links.length) {
    console.log("link length diff");
    console.log(publishedCollection.links.length);
    console.log(revision.links.length);
    return true;
  }
  //Check links for differences
  if (checkListForChanges(revision.links, publishedCollection.links)) {
    console.log("changes in link lists");
    return true;
  }
  if (
    checkDatasetsForChanges(revision.datasets, publishedCollection.datasets)
  ) {
    console.log("datasets changed");
    return true;
  }
  console.log("no changes!");
  return false;
}
