import {
  useMutation,
  useQuery,
  useQueryClient,
  UseQueryResult,
} from "react-query";
import { Collection, VISIBILITY_TYPE } from "src/common/entities";
import { buildSummaryCitation } from "src/common/queries/filter";
import { apiTemplateToUrl } from "src/common/utils/apiTemplateToUrl";
import { API_URL } from "src/configs/configs";
import { API } from "../API";
import HTTP_STATUS_CODE from "../constants/HTTP_STATUS_CODE";
import checkForRevisionChange from "../utils/checkForRevisionChange";
import { isTombstonedCollection } from "../utils/typeGuards";
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

export function useCollections(): UseQueryResult<CollectionResponsesMap> {
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

function generateDatasetMap(json: any) {
  const datasetMap = new Map() as Collection["datasets"];
  for (const dataset of json.datasets) {
    datasetMap.set(dataset.original_id || dataset.id, dataset);
  }
  return datasetMap;
}

export interface TombstonedCollection {
  tombstone: true;
}

function fetchCollection(allCollections: CollectionResponsesMap | undefined) {
  return async function (
    id: string,
    visibility: VISIBILITY_TYPE
  ): Promise<Collection | TombstonedCollection | null> {
    if (!id) {
      return null;
    }

    const baseUrl = apiTemplateToUrl(API_URL + API.COLLECTION, { id });

    const finalUrl =
      visibility === VISIBILITY_TYPE.PRIVATE
        ? baseUrl + `?visibility=${VISIBILITY_TYPE.PRIVATE}`
        : baseUrl;

    let response = await fetch(finalUrl, DEFAULT_FETCH_OPTIONS);
    let json = await response.json();

    if (response.status === HTTP_STATUS_CODE.GONE) {
      return { tombstone: true };
    }
    if (!response.ok) {
      throw json;
    }

    const collection = { ...json, datasets: generateDatasetMap(json) };

    let publishedCounterpart;

    if (allCollections) {
      const collectionsWithID = allCollections.get(id);

      collection.has_revision = collectionsWithID && collectionsWithID.size > 1;
    }

    if (visibility === VISIBILITY_TYPE.PRIVATE && collection.has_revision) {
      response = await fetch(baseUrl, DEFAULT_FETCH_OPTIONS);
      json = await response.json();

      if (response.ok) {
        publishedCounterpart = {
          ...json,
          datasets: generateDatasetMap(json),
        };
      }
    }

    // check for diffs between revision and published collection
    if (collection.has_revision && visibility === VISIBILITY_TYPE.PRIVATE) {
      collection.revision_diff = checkForRevisionChange(
        collection,
        publishedCounterpart
      );
    }

    // Add summary citation to collection.
    collection.summaryCitation = buildSummaryCitation(
      collection.publisher_metadata
    );

    return collection;
  };
}

export function useCollection({
  id = "",
  visibility = VISIBILITY_TYPE.PUBLIC,
}: {
  id?: string;
  visibility: VISIBILITY_TYPE;
}): UseQueryResult<Collection | TombstonedCollection | null> {
  const { data: collections } = useCollections();
  const queryFn = fetchCollection(collections);

  return useQuery<Collection | TombstonedCollection | null>(
    [USE_COLLECTION, id, visibility, collections],
    () => queryFn(id, visibility),
    { enabled: !!collections }
  );
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
  const queryClient = useQueryClient();

  return useMutation(createCollection, {
    onSuccess: () => {
      queryClient.invalidateQueries([USE_COLLECTIONS]);
    },
  });
}

export function formDataToObject(formData: FormData): {
  [key: string]: unknown;
} {
  const payload: { [key: string]: unknown } = {};

  formData.forEach((value, key: string) => {
    const translatedKey = key.replace("-", "_");

    payload[translatedKey] = value;
  });

  return payload;
}

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
  const queryCache = useQueryClient();

  return useMutation(collectionUploadLinks, {
    onSuccess: () => {
      queryCache.invalidateQueries([USE_COLLECTION, id, visibility]);
    },
  });
}

async function deleteCollection({
  collectionID,
  visibility = VISIBILITY_TYPE.PRIVATE,
}: {
  collectionID: Collection["id"];
  visibility?: VISIBILITY_TYPE;
}) {
  const baseUrl = apiTemplateToUrl(API_URL + API.COLLECTION, {
    id: collectionID,
  });
  const finalUrl = baseUrl + `?visibility=${visibility}`;

  const response = await fetch(finalUrl, DELETE_FETCH_OPTIONS);

  if (!response.ok) {
    throw await response.json();
  }
}

export function useDeleteCollection(
  id = "",
  visibility = VISIBILITY_TYPE.PRIVATE
) {
  const queryClient = useQueryClient();

  return useMutation(deleteCollection, {
    onError: (
      _,
      __,
      context: { previousCollections: CollectionResponsesMap } | undefined
    ) => {
      queryClient.setQueryData([USE_COLLECTIONS], context?.previousCollections);
    },
    onMutate: async () => {
      await queryClient.cancelQueries([USE_COLLECTIONS]);

      const previousCollections = queryClient.getQueryData([
        USE_COLLECTIONS,
      ]) as CollectionResponsesMap;

      const newCollections = new Map(previousCollections);
      const collectionsWithID = newCollections.get(id);
      // If we're deleting a public collection or there is no revision, nuke it from the cache
      if (
        visibility === VISIBILITY_TYPE.PUBLIC ||
        (collectionsWithID && collectionsWithID.entries.length > 1)
      ) {
        newCollections.delete(id);
      } else {
        // Otherwise, we need to preserve the public collection
        collectionsWithID?.delete(VISIBILITY_TYPE.PRIVATE);
        if (collectionsWithID) newCollections.set(id, collectionsWithID);
      }
      queryClient.setQueryData([USE_COLLECTIONS], newCollections);

      return { previousCollections };
    },
    onSuccess: () => {
      return Promise.all([
        queryClient.invalidateQueries([USE_COLLECTIONS]),
        queryClient.removeQueries([USE_COLLECTION, id, visibility], {
          exact: false,
        }),
      ]);
    },
  });
}

export type PublishCollection = {
  id: Collection["id"];
  payload: string;
};

async function publishCollection({ id, payload }: PublishCollection) {
  const url = apiTemplateToUrl(API_URL + API.COLLECTION_PUBLISH, { id });

  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    body: payload,
    method: "POST",
  });

  if (!response.ok) {
    throw await response.json();
  }
}

export function usePublishCollection() {
  const queryClient = useQueryClient();

  return useMutation(publishCollection, {
    onSuccess: () => {
      // (thuang): We don't need to invalidate `[USE_COLLECTION, id, visibility]`
      // because `visibility` has changed from PRIVATE to PUBLIC
      queryClient.invalidateQueries([USE_COLLECTIONS]);
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

export function useEditCollection(collectionID?: Collection["id"]) {
  const queryClient = useQueryClient();

  const { data: collections } = useCollections();

  const { data: collection } = useCollection({
    id: collectionID,
    visibility: VISIBILITY_TYPE.PRIVATE,
  });

  const { data: publishedCollection } = useCollection({
    id: collectionID,
    visibility: VISIBILITY_TYPE.PUBLIC,
  });

  return useMutation(editCollection, {
    onSuccess: (newCollection) => {
      queryClient.setQueryData(
        [USE_COLLECTION, collectionID, VISIBILITY_TYPE.PRIVATE, collections],
        () => {
          let revision_diff;
          if (isTombstonedCollection(newCollection)) {
            return newCollection;
          } else if (
            !isTombstonedCollection(collection) &&
            !isTombstonedCollection(publishedCollection) &&
            collection?.has_revision &&
            publishedCollection
          ) {
            revision_diff = checkForRevisionChange(
              newCollection,
              publishedCollection
            );

            return { ...collection, ...newCollection, revision_diff };
          }
        }
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
  return response.json();
};

export function useCreateRevision(callback: () => void) {
  const queryClient = useQueryClient();

  return useMutation(createRevision, {
    onSuccess: async (collection: Collection) => {
      await queryClient.invalidateQueries([USE_COLLECTIONS]);
      await queryClient.invalidateQueries([USE_COLLECTION, collection.id]);
      callback();
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
  const queryClient = useQueryClient();

  return useMutation(reuploadDataset, {
    onSuccess: () => {
      queryClient.invalidateQueries([
        USE_COLLECTION,
        collectionId,
        VISIBILITY_TYPE.PRIVATE,
      ]);
    },
  });
}
