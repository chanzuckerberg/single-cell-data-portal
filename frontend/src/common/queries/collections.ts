import { useMutation, useQuery, useQueryCache } from "react-query";
import { Collection, Dataset, VISIBILITY_TYPE } from "src/common/entities";
import { apiTemplateToUrl } from "src/common/utils/apiTemplateToUrl";
import { API_URL } from "src/configs/configs";
import { API } from "../API";
import {
  DEFAULT_FETCH_OPTIONS,
  DELETE_FETCH_OPTIONS,
  JSON_BODY_FETCH_OPTIONS,
} from "./common";
import { ENTITIES } from "./entities";

export const USE_COLLECTIONS = {
  entities: [ENTITIES.COLLECTION],
  id: "collections",
};

function idError(id: string | null) {
  if (!id) {
    throw Error("No collection id given");
  }
}

export enum REVISION_STATUS {
  NOT_STARTED,
  STARTED,
  DISABLED,
}

export interface CollectionResponse {
  id: string;
  created_at: number;
  visibility: VISIBILITY_TYPE;
}

export interface RevisionResponse extends CollectionResponse {
  revision: REVISION_STATUS;
}

// if we return as a Map it will be easier to fetch collections by id
export type CollectionResponsesMap = Map<
  CollectionResponse["id"],
  Map<VISIBILITY_TYPE, CollectionResponse>
>;

async function fetchCollections(): Promise<CollectionResponsesMap> {
  const json = await (
    await fetch(API_URL + API.COLLECTIONS, DEFAULT_FETCH_OPTIONS)
  ).json();

  const collectionsMap: CollectionResponsesMap = new Map();

  for (const collection of json.collections as CollectionResponse[]) {
    const collectionsWithId = collectionsMap.get(collection.id) || new Map();
    collectionsWithId.set(collection.visibility, collection);
    collectionsMap.set(collection.id, collectionsWithId);
  }

  return collectionsMap;
}

export function useCollections() {
  return useQuery([USE_COLLECTIONS], fetchCollections);
}

export const USE_COLLECTION = {
  entities: [ENTITIES.COLLECTION, ENTITIES.DATASET],
  id: "collection",
};

export type CollectionError = {
  detail: string;
  status: number;
  title: string;
  type: string;
};

async function fetchCollection(
  _: unknown,
  id: string,
  visibility: VISIBILITY_TYPE
): Promise<Collection | null> {
  if (!id) {
    return null;
  }

  const baseUrl = apiTemplateToUrl(API_URL + API.COLLECTION, { id });

  const finalUrl =
    visibility === VISIBILITY_TYPE.PRIVATE
      ? baseUrl + `?visibility=${VISIBILITY_TYPE.PRIVATE}`
      : baseUrl;

  const response = await fetch(finalUrl, DEFAULT_FETCH_OPTIONS);
  const result = await response.json();

  if (!response.ok) {
    throw result;
  }

  return result;
}

const ignoredCollectionFields = [
  "visibility",
  "created_at",
  "updated_at",
  "is_revision",
  "revision_diff",
] as Array<keyof Collection>;
const ignoredDatasetFields = [
  "created_at",
  "updated_at",
  "collection_visibility",
  "original_uuid",
] as Array<keyof Dataset>;

function useCollectionFetch({ id = "", visibility = VISIBILITY_TYPE.PUBLIC }) {
  const { data: collections } = useCollections();
  const queryCache = useQueryCache();

  return useQuery<Collection | null>(
    [USE_COLLECTION, id, visibility, "FETCH"],
    fetchCollection,
    {
      onSuccess: (data) => {
        // Check if there are multiple collections with this ID
        if (collections && data) {
          // Set is_revision to true or false
          queryCache.setQueryData([USE_COLLECTION, id, visibility], () => {
            const newData = data;
            const collectionsWithID = collections.get(data.id);
            if (!collectionsWithID) return newData;
            newData.is_revision = collectionsWithID.size > 1;
            if (newData.is_revision && visibility === VISIBILITY_TYPE.PRIVATE) {
              const publishedCollection = queryCache.getQueryData([
                USE_COLLECTION,
                id,
                visibility,
                VISIBILITY_TYPE.PUBLIC,
              ]) as Collection;
              if (!publishCollection) return newData;
              let revisionChange = false;

              // For some reason key would be typed as "string" if I didn't explicitly type it here
              let collectionKey = "" as keyof Collection;
              for (collectionKey in publishedCollection) {
                if (
                  !revisionChange &&
                  !ignoredCollectionFields.includes(collectionKey)
                ) {
                  revisionChange =
                    publishedCollection[collectionKey] !==
                    newData[collectionKey];
                }
              }
              publishedCollection.datasets.forEach((publishedDataset) => {
                if (revisionChange) return newData;
                let datasetKey = "" as keyof Dataset;
                for (datasetKey in publishedDataset) {
                  if (
                    revisionChange &&
                    ignoredDatasetFields.includes(datasetKey)
                  ) {
                    revisionChange =
                      publishedDataset[datasetKey] !== newDataset[datasetKey];
                  }
                }
              });
            }

            return newData;
          });
        }
      },
    }
  );
}

export function useCollection({
  id = "",
  visibility = VISIBILITY_TYPE.PUBLIC,
}) {
  useCollectionFetch({ id, visibility });
  return useQuery<Collection | null>([USE_COLLECTION, id, visibility, "FETCH"]);
}

export async function createCollection(payload: string): Promise<string> {
  const response = await fetch(`${API_URL}${API.CREATE_COLLECTION}`, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    body: payload,

    method: "POST",
  });

  const json = await response.json();

  if (!response.ok) {
    throw json;
  }

  return json.collection_uuid;
}

export function useCreateCollection() {
  const queryCache = useQueryCache();

  return useMutation(createCollection, {
    onSuccess: () => {
      queryCache.invalidateQueries(USE_COLLECTIONS);
    },
  });
}

export const formDataToObject = function (formData: FormData) {
  const payload: { [key: string]: unknown } = {};

  formData.forEach((value, key: string) => {
    const translatedKey = key.replace("-", "_");

    payload[translatedKey] = value;
  });

  return payload;
};

export type CollectionUploadLinks = {
  collectionId: string;
  payload: string;
};

async function collectionUploadLinks({
  collectionId,
  payload,
}: CollectionUploadLinks) {
  const url = apiTemplateToUrl(API_URL + API.COLLECTION_UPLOAD_LINKS, {
    id: collectionId,
  });

  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    body: payload,
    method: "POST",
  });

  const json = await response.json();

  if (!response.ok) {
    throw json;
  }

  return json.dataset_uuid;
}

export function useCollectionUploadLinks(
  id: string,
  visibility: VISIBILITY_TYPE
) {
  const queryCache = useQueryCache();

  return useMutation(collectionUploadLinks, {
    onSuccess: () => {
      queryCache.invalidateQueries([USE_COLLECTION, id, visibility]);
    },
  });
}

async function deleteCollection(collectionID: Collection["id"]) {
  const baseUrl = apiTemplateToUrl(API_URL + API.COLLECTION, {
    id: collectionID,
  });
  const finalUrl = baseUrl + `?visibility=${VISIBILITY_TYPE.PRIVATE}`;

  const response = await fetch(finalUrl, DELETE_FETCH_OPTIONS);

  if (!response.ok) {
    throw await response.json();
  }
}

export function useDeleteCollection(id = "") {
  if (!id) {
    throw new Error("No collection id given");
  }

  const queryCache = useQueryCache();

  return useMutation(deleteCollection, {
    onSuccess: () => {
      queryCache.invalidateQueries([USE_COLLECTIONS]);
    },
  });
}

async function publishCollection(id: Collection["id"]) {
  const url = apiTemplateToUrl(API_URL + API.COLLECTION_PUBLISH, { id });

  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    method: "POST",
  });

  if (!response.ok) {
    throw await response.json();
  }
}

export function usePublishCollection() {
  const queryCache = useQueryCache();

  return useMutation(publishCollection, {
    onSuccess: () => {
      // (thuang): We don't need to invalidate `[USE_COLLECTION, id, visibility]`
      // because `visibility` has changed from PRIVATE to PUBLIC
      queryCache.invalidateQueries([USE_COLLECTIONS]);
    },
  });
}

const editCollection = async function ({
  id,
  payload,
}: {
  id: string;
  payload: string;
}): Promise<Collection> {
  idError(id);
  if (!payload) {
    throw Error("No payload given");
  }

  const url = apiTemplateToUrl(API_URL + API.COLLECTION, { id });

  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    body: payload,
    method: "PUT",
  });

  const result = await response.json();

  if (!response.ok) throw result;
  return result;
};

export function useEditCollection() {
  const queryCache = useQueryCache();

  return useMutation(editCollection, {
    onSuccess: (collection: Collection) => {
      return queryCache.setQueryData(
        [USE_COLLECTION, collection.id, collection.visibility],
        collection
      );
    },
  });
}

const createRevision = async function (id: string) {
  idError(id);
  const url = apiTemplateToUrl(API_URL + API.COLLECTION, { id });

  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    method: "POST",
  });

  if (!response.ok) throw await response.json();
};

export function useCreateRevision(callback: () => void) {
  const queryCache = useQueryCache();

  return useMutation(createRevision, {
    onSuccess: () => {
      callback();
      return queryCache.invalidateQueries([USE_COLLECTIONS]);
    },
  });
}

export interface ReuploadLink {
  payload: string;
  collectionId: string;
}

const reuploadDataset = async function ({
  payload,
  collectionId,
}: ReuploadLink) {
  if (!payload) {
    throw Error("No payload given");
  }
  idError(collectionId);

  const url = apiTemplateToUrl(API_URL + API.COLLECTION_UPLOAD_LINKS, {
    id: collectionId,
  });

  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    body: payload,
    method: "PUT",
  });

  const result = await response.json();
  if (!response.ok) throw result;

  return result.dataset_uuid;
};

export function useReuploadDataset(collectionId: string) {
  const queryCache = useQueryCache();

  return useMutation(reuploadDataset, {
    onSuccess: () => {
      queryCache.invalidateQueries([
        USE_COLLECTION,
        collectionId,
        VISIBILITY_TYPE.PRIVATE,
      ]);
    },
  });
}
