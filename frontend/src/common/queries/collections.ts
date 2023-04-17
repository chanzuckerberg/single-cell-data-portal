import {
  useMutation,
  UseMutationResult,
  useQuery,
  useQueryClient,
  UseQueryResult,
} from "react-query";
import {
  Collection,
  COLLECTION_STATUS,
  VISIBILITY_TYPE,
} from "src/common/entities";
import {
  buildSummaryCitation,
  ProcessedCollectionResponse,
  USE_COLLECTIONS_INDEX,
  USE_DATASETS_INDEX,
} from "src/common/queries/filter";
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
import { QueryClient } from "react-query/core";

/**
 * Never expire cached collection / collections.
 */
const DEFAULT_QUERY_OPTIONS = {
  staleTime: Infinity,
};

/**
 * Error text returned from BE when DOI format is identified as invalid.
 */
const INVALID_DOI_FORMAT_MESSAGE = "Invalid DOI";

/**
 * Error text returned from BE when DOI is identified as invalid, that is, DOI can not be found on Crossref.
 */
const INVALID_DOI_MESSAGE = "DOI cannot be found on Crossref";

/**
 * Model returned from create collection: collection ID if create was successful, otherwise flag indicating invalid DOI.
 */
export interface CollectionCreateResponse {
  collectionId?: string;
  isInvalidDOI?: boolean;
}

/**
 * Model returned from edit collection: the collection itself if edit was successful, otherwise flag indicating invalid
 * DOI.
 */
export interface CollectionEditResponse {
  collection?: Collection;
  isInvalidDOI?: boolean;
}

/**
 * Possible set of error keys returned from the BE.
 */
enum ERROR_KEY {
  "DOI" = "link_type",
}

/**
 * Possible set of error values returned from the BE.
 */
enum ERROR_VALUE {
  "DOI" = "DOI",
}

/**
 * Error information returned from create and edit collection API endpoints. For example,
 * see the single element in the "detail" array in the following error response:
 *
 * {
 *   "detail": [{
 *       link_type: "DOI",
 *       reason: "Invalid DOI"
 *   }],
 *   "status": 400,
 *   "title": "Bad Request",
 *   "type": "about:blank"
 * }
 */
type Error = { [key in ERROR_KEY]?: ERROR_VALUE } & { reason: string };

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
  revision_of?: string;
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
    const publicID = collection.revision_of || collection.id;
    const collectionsWithID = collectionsMap.get(publicID) || new Map();
    collectionsWithID.set(collection.visibility, collection);
    collectionsMap.set(publicID, collectionsWithID);
  }

  return collectionsMap;
}

export function useCollections(): UseQueryResult<CollectionResponsesMap> {
  return useQuery([USE_COLLECTIONS], fetchCollections, {
    ...DEFAULT_QUERY_OPTIONS,
  });
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

function fetchCollection() {
  return async function (
    id: string
  ): Promise<Collection | TombstonedCollection | null> {
    if (!id) {
      return null;
    }

    const url = apiTemplateToUrl(API_URL + API.COLLECTION, { id });

    let response = await fetch(url, DEFAULT_FETCH_OPTIONS);
    let json = await response.json();

    if (response.status === HTTP_STATUS_CODE.GONE) {
      return { tombstone: true };
    }
    if (!response.ok) {
      throw json;
    }

    const collection: Collection = {
      ...json,
      datasets: generateDatasetMap(json),
    };

    let publishedCounterpart;

    if (collection.revision_of) {
      const publicCollectionURL = apiTemplateToUrl(API_URL + API.COLLECTION, {
        id: collection.revision_of,
      });
      response = await fetch(publicCollectionURL, DEFAULT_FETCH_OPTIONS);
      json = await response.json();

      if (response.ok) {
        publishedCounterpart = {
          ...json,
          datasets: generateDatasetMap(json),
        };
      }
    }

    // check for diffs between revision and published collection
    if (collection.revision_of) {
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
}: {
  id?: string;
}): UseQueryResult<Collection | TombstonedCollection | null> {
  const queryFn = fetchCollection();
  return useQuery<Collection | TombstonedCollection | null>(
    [USE_COLLECTION, id],
    () => queryFn(id),
    {
      ...DEFAULT_QUERY_OPTIONS,
    }
  );
}

export async function createCollection(
  payload: string
): Promise<CollectionCreateResponse> {
  const response = await fetch(`${API_URL}${API.CREATE_COLLECTION}`, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    body: payload,

    method: "POST",
  });

  const json = await response.json();

  // Check for validation errors. Currently only DOI is validated by the BE; this can be generalized once all fields
  // are validated by the BE.
  if (isInvalidDOI(response.status, json.detail)) {
    return { isInvalidDOI: true };
  }

  if (!response.ok) {
    throw json;
  }

  return {
    collectionId: json.collection_id,
  };
}

export function useCreateCollection() {
  const queryClient = useQueryClient();
  return useMutation(createCollection, {
    onSuccess: async () => {
      await queryClient.invalidateQueries([USE_COLLECTIONS_INDEX]);
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

  return json.dataset_id;
}

export function useCollectionUploadLinks(id: string) {
  const queryCache = useQueryClient();

  return useMutation(collectionUploadLinks, {
    onSuccess: () => {
      queryCache.invalidateQueries([USE_COLLECTION, id]);
    },
  });
}

async function deleteCollection({
  collectionID,
}: {
  collectionID: Collection["id"];
}) {
  const finalURL = apiTemplateToUrl(API_URL + API.COLLECTION, {
    id: collectionID,
  });

  const response = await fetch(finalURL, DELETE_FETCH_OPTIONS);

  if (!response.ok) {
    throw await response.json();
  }
}

export function useDeleteCollection(
  id = "",
  visibility = ""
): UseMutationResult<
  void,
  unknown,
  { collectionID: string },
  { previousCollections: CollectionResponsesMap }
> {
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
        queryClient.removeQueries([USE_COLLECTION, id], {
          exact: false,
        }),
      ]);
    },
  });
}

export type PublishCollection = {
  collection: Collection;
  payload: string;
};

async function publishCollection({
  collection,
  payload,
}: PublishCollection): Promise<Collection> {
  const url = apiTemplateToUrl(API_URL + API.COLLECTION_PUBLISH, {
    id: collection.id,
  });
  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    body: payload,
    method: "POST",
  });
  if (!response.ok) {
    throw await response.json();
  }
  return collection;
}

export function usePublishCollection() {
  const queryClient = useQueryClient();
  return useMutation(publishCollection, {
    onSuccess: async (collection: Collection): Promise<void> => {
      await queryClient.invalidateQueries([USE_COLLECTIONS_INDEX]);
      await queryClient.prefetchQuery([USE_COLLECTIONS_INDEX]);
      await queryClient.invalidateQueries([USE_DATASETS_INDEX]);
      await queryClient.prefetchQuery([USE_DATASETS_INDEX]);
      await queryClient.invalidateQueries([USE_COLLECTION, collection.id]);
      if (collection.revision_of) {
        await queryClient.invalidateQueries([
          USE_COLLECTION,
          collection.revision_of,
        ]);
      }
    },
  });
}

const editCollection = async function ({
  id,
  payload,
}: {
  id: string;
  payload: string;
}): Promise<CollectionEditResponse> {
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

  // Check for validation errors. Currently only DOI is validated by the BE; this can be generalized once all fields
  // are validated by the BE.
  if (isInvalidDOI(response.status, result.detail)) {
    return { isInvalidDOI: true };
  }

  if (!response.ok) {
    throw result;
  }

  return {
    collection: result,
  };
};

export function useEditCollection(
  collectionID?: Collection["id"],
  publicID?: Collection["id"]
): UseMutationResult<
  CollectionEditResponse,
  unknown,
  { id: string; payload: string },
  unknown
> {
  const queryClient = useQueryClient();
  const { data: collection } = useCollection({
    //revision
    id: collectionID,
  });
  const { data: publishedCollection } = useCollection({
    //published collection
    id: publicID,
  });

  return useMutation(editCollection, {
    // newCollection is the result of the PUT on the revision
    onSuccess: async ({ collection: newCollection }) => {
      // Check for updated collection: it's possible server-side validation errors have occurred where the error has
      // been swallowed (allowing error messages to be displayed on the edit form) and success flow is executed even
      // though update did not occur.
      if (!newCollection) {
        return;
      }
      queryClient.setQueryData([USE_COLLECTION, collectionID], () => {
        let revision_diff;
        if (isTombstonedCollection(newCollection)) {
          return newCollection;
        } else if (
          !isTombstonedCollection(collection) &&
          !isTombstonedCollection(publishedCollection) &&
          collection?.revision_of &&
          publishedCollection
        ) {
          revision_diff = checkForRevisionChange(
            newCollection,
            publishedCollection
          );

          return { ...collection, ...newCollection, revision_diff };
        }
        return { ...collection, ...newCollection };
      });
      await queryClient.invalidateQueries([USE_COLLECTIONS_INDEX]);
      await queryClient.prefetchQuery([USE_COLLECTIONS_INDEX]);
      await queryClient.invalidateQueries([USE_DATASETS_INDEX]);
      await queryClient.prefetchQuery([USE_DATASETS_INDEX]);
    },
  });
}

const createRevision = async function (id: Collection["id"]) {
  idError(id);
  const url = apiTemplateToUrl(API_URL + API.COLLECTION, { id });

  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    method: "POST",
  });

  if (!response.ok) throw await response.json();
  return response.json();
};

/**
 * Sets USE_COLLECTION, updates the public collection with:
 * - revision_diff: the revision difference between the public collection and the private revision, and
 * - revising_in: the id of the private revision collection.
 * @param revision - Private revision collection.
 * @param queryClient - QueryClient to update the cache.
 */
function createRevisionSetUseCollection(
  revision: Collection,
  queryClient: QueryClient
): void {
  // Grab the public collection cache.
  const publicCollection = queryClient.getQueryData([
    USE_COLLECTION,
    revision.revision_of,
  ]) as Collection;
  // Update the public collection cache.
  queryClient.setQueryData([USE_COLLECTION, revision.revision_of], {
    ...publicCollection,
    revision_diff: false,
    revising_in: revision.id,
  });
  // Create the private revision collection cache.
  queryClient.setQueryData([USE_COLLECTION, revision.id], revision);
}

/**
 * Sets USE_COLLECTIONS_INDEX, updates the public collection of the new private revision with:
 * - revisedBy: the id of the new private revision, and
 * - status: the public collection status.
 * @param revision - Private revision collection.
 * @param queryClient - QueryClient to update the cache.
 */
function createRevisionSetUseCollectionsIndex(
  revision: Collection,
  queryClient: QueryClient
): void {
  if (!revision.revision_of) {
    throw Error("Create revision response does not have a revision_of field");
  }
  const updatedCollectionsById = new Map(
    queryClient.getQueryData([USE_COLLECTIONS_INDEX])
  ) as Map<string, ProcessedCollectionResponse>;
  const collection = updatedCollectionsById.get(revision.revision_of);
  const updatedCollection = {
    ...collection,
    revisedBy: revision.id,
    status: [COLLECTION_STATUS.PUBLISHED, COLLECTION_STATUS.REVISION],
  } as ProcessedCollectionResponse;
  updatedCollectionsById.set(revision.revision_of, updatedCollection);
  queryClient.setQueryData([USE_COLLECTIONS_INDEX], updatedCollectionsById);
}

export function useCreateRevision() {
  const queryClient = useQueryClient();
  return useMutation(createRevision, {
    onSuccess: async (revision) => {
      createRevisionSetUseCollectionsIndex(revision, queryClient);
      createRevisionSetUseCollection(revision, queryClient);
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

  return result.dataset_id;
};

export function useReuploadDataset(
  collectionId: string
): UseMutationResult<unknown, unknown, ReuploadLink, unknown> {
  const queryClient = useQueryClient();

  return useMutation(reuploadDataset, {
    onSuccess: () => {
      queryClient.invalidateQueries([USE_COLLECTION, collectionId]);
    },
  });
}

/**
 * Determine if a submitted DOI has failed validation on the BE.
 *
 * Expected response for invalid DOI:
 * {
 *   "detail": [{
 *       link_type: "DOI",
 *       reason: "DOI cannot be found on Crossref"
 *   }],
 *   "status": 400,
 *   "title": "Bad Request",
 *   "type": "about:blank"
 * }
 *
 * Expected response for DOI with an invalid format:
 * {
 *   "detail": [{
 *       link_type: "DOI",
 *       reason: "Invalid DOI"
 *   }],
 *   "status": 400,
 *   "title": "Bad Request",
 *   "type": "about:blank"
 * }
 *
 * TODO generalize beyond DOI link type once all links are validated on the BE (#1916).
 *
 * @param status - Response status returned from server.
 * @param errors - Array of errors returned from server, if any.
 * @returns True if DOI has been identified as invalid by the BE.
 */
function isInvalidDOI(status: number, errors?: Error[]): boolean {
  if (status !== HTTP_STATUS_CODE.BAD_REQUEST || !errors) {
    return false;
  }

  // Check if the errors returned from the server contain a DOI error.
  const doiError = findErrorByKey(errors, ERROR_KEY.DOI);
  if (!doiError) {
    return false;
  }

  // There's a DOI error; check if it's an error we report on.
  const { reason } = doiError;
  return (
    reason === INVALID_DOI_MESSAGE || reason === INVALID_DOI_FORMAT_MESSAGE
  );
}

/**
 * Find the error with the given key.
 * @param errors - Array of errors returned from server.
 * @param key - Key of error to find.
 * @returns The error with the given key, if any.
 */
function findErrorByKey(errors: Error[], key: ERROR_KEY): Error | undefined {
  return errors.find((error) => {
    return Object.keys(error).find((errorKey) => errorKey === key);
  });
}
