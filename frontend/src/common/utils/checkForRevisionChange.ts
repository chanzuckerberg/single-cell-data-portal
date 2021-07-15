import { isEqual } from "lodash";
import { Collection, Dataset } from "../entities";

const IGNORED_COLLECTION_FIELDS = [
  "visibility",
  "created_at",
  "updated_at",
  "is_revision",
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
  "processing_status",
  "dataset_assets",
] as Array<keyof Dataset>;

function checkListForChanges(
  revisedList: Array<unknown>,
  publishedList: Array<unknown>
) {
  if (!revisedList || publishedList.length !== revisedList.length) {
    return true;
  }
  if (publishedList.length > 0 && publishedList[0] instanceof Object) {
    for (let index = 0; index < publishedList.length; index++) {
      const publishedObj = publishedList[index] as Record<string, unknown>;
      if (!isEqual(publishedObj, revisedList[index])) {
        return true;
      }
    }
  }
}

function checkDatasetsForChanges(
  revisedDatasets: Map<Dataset["id"], Dataset>,
  publishedDatasets: Map<Dataset["id"], Dataset>
) {
  if (publishedDatasets.size !== revisedDatasets.size) return true;
  // Check dataset fields for differences
  Array.from(publishedDatasets.values()).forEach((publishedDataset) => {
    const revisedDataset =
      revisedDatasets.get(publishedDataset.id) || ({} as Dataset);
    let datasetKey = "" as keyof Dataset;
    for (datasetKey in publishedDataset) {
      if (IGNORED_DATASET_FIELDS.includes(datasetKey)) {
        continue; // ignore these fields (I hate using continue, but sonarjs considers this best)
      } else if (publishedDataset[datasetKey] instanceof Array) {
        // Dataset entry is an array
        if (
          checkListForChanges(
            publishedDataset[datasetKey] as Array<unknown>,
            revisedDataset[datasetKey] as Array<unknown>
          )
        ) {
          return true;
        }
      } else if (
        // Dataset entry is an object
        publishedDataset[datasetKey] instanceof Object &&
        !isEqual(publishedDataset[datasetKey], revisedDataset[datasetKey])
      ) {
        return true;
      } else if (
        // scalar value
        publishedDataset[datasetKey] !== revisedDataset[datasetKey]
      ) {
        return true;
      }
    }
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
      return true;
    }
  }
  if (publishedCollection.links.length !== revision.links.length) return true;
  //Check links for differences
  publishedCollection.links.forEach((link, index) => {
    if (link !== revision.links[index]) {
      return true;
    }
  });

  if (checkDatasetsForChanges(revision.datasets, publishedCollection.datasets))
    return true;

  return false;
}
