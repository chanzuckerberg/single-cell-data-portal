import {
  useMutation,
  UseMutationResult,
  useQuery,
  useQueryClient,
  UseQueryResult,
} from "react-query";
import { Collection, COLLECTION_STATUS } from "src/common/entities";
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
 * Never expire cached query.
 */
const DEFAULT_QUERY_OPTIONS = {
  staleTime: Infinity,
};

/**
 * Cached query matching the refetch predicate, that are not being rendered, will be invalidated and refetched
 * in the background.
 */
const DEFAULT_BACKGROUND_REFETCH = {
  refetchInactive: true,
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

function idError(id: string | null) {
  if (!id) {
    throw Error("No collection id given");
  }
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
    onSuccess: async (): Promise<void> => {
      await queryClient.invalidateQueries(
        [USE_COLLECTIONS_INDEX],
        DEFAULT_BACKGROUND_REFETCH
      );
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
    onSuccess: async (): Promise<void> => {
      await queryCache.invalidateQueries([USE_COLLECTION, id]);
    },
  });
}

async function deleteCollection(collection: Collection): Promise<Collection> {
  const finalURL = apiTemplateToUrl(API_URL + API.COLLECTION, {
    id: collection.id,
  });

  const response = await fetch(finalURL, DELETE_FETCH_OPTIONS);

  if (!response.ok) {
    throw await response.json();
  }
  return collection;
}

export function useDeleteCollection(): UseMutationResult<
  Collection,
  unknown,
  Collection
> {
  const queryClient = useQueryClient();
  return useMutation(deleteCollection, {
    onSuccess: async (collection: Collection): Promise<void> => {
      queryClient.removeQueries([USE_COLLECTION, collection.id]);
      await queryClient.invalidateQueries(
        [USE_COLLECTIONS_INDEX],
        DEFAULT_BACKGROUND_REFETCH
      );
      await queryClient.invalidateQueries(
        [USE_DATASETS_INDEX],
        DEFAULT_BACKGROUND_REFETCH
      );
      if (collection.revision_of) {
        await queryClient.invalidateQueries(
          [USE_COLLECTION, collection.revision_of],
          DEFAULT_BACKGROUND_REFETCH
        );
      }
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
      await queryClient.invalidateQueries(
        [USE_COLLECTIONS_INDEX],
        DEFAULT_BACKGROUND_REFETCH
      );
      await queryClient.invalidateQueries(
        [USE_DATASETS_INDEX],
        DEFAULT_BACKGROUND_REFETCH
      );
      // Invalidate the private or private revision collection.
      // If the collection is a private revision, the query should be marked as invalid without executing a background
      // re-fetch. A background refresh on the "active" private revision collection query will trigger an immediate
      // re-fetch of the collection. However, the response returns published values where the "id" corresponds
      // to the published collection's "id", which is the "revision_of" value of the private revision.
      // Additionally, the visibility of the collection is updated to "PUBLIC".
      // This re-fetch happens before the "usePublishCollection" mutate function executes the "onSuccess" callback,
      // resulting in unintended consequences. For example, the "Publish" button unmounts before the user is routed
      // to the updated collection because the visibility is no longer "PRIVATE."
      const inRevision = Boolean(collection.revision_of);
      await queryClient.invalidateQueries([USE_COLLECTION, collection.id], {
        refetchActive: !inRevision, // If the collection is in revision, mark as invalid without executing a re-fetch.
      });
      if (inRevision) {
        // If the collection is in revision, invalidate the published collection.
        await queryClient.invalidateQueries(
          [USE_COLLECTION, collection.revision_of],
          DEFAULT_BACKGROUND_REFETCH
        );
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
    onSuccess: async ({ collection: newCollection }): Promise<void> => {
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
      await queryClient.invalidateQueries(
        [USE_COLLECTIONS_INDEX],
        DEFAULT_BACKGROUND_REFETCH
      );
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
    onSuccess: async (revision): Promise<void> => {
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
    onSuccess: async () => {
      await queryClient.invalidateQueries([USE_COLLECTION, collectionId]);
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
