import { useMemo } from "react";
import { QueryStatus, useQuery, UseQueryResult } from "react-query";
import { API } from "src/common/API";
import {
  ACCESS_TYPE,
  Author,
  Collection,
  COLLECTION_STATUS,
  Consortium,
  IS_PRIMARY_DATA,
  Ontology,
  PublisherMetadata,
  VISIBILITY_TYPE,
} from "src/common/entities";
import {
  buildExplicitOntologyTermId,
  buildInferredOntologyTermId,
} from "src/common/hooks/useCategoryFilter/common/multiPanelOntologyUtils";
import { DEFAULT_FETCH_OPTIONS } from "src/common/queries/common";
import { ENTITIES } from "src/common/queries/entities";
import {
  COLLATOR_CASE_INSENSITIVE,
  PUBLICATION_DATE_VALUES,
  SELF_REPORTED_ETHNICITY_DENY_LIST,
  SUSPENSION_TYPE_DENY_LIST,
} from "src/components/common/Filter/common/constants";
import {
  Categories,
  CollectionRow,
  DatasetRow,
} from "src/components/common/Filter/common/entities";
import { checkIsOverMaxCellCount } from "src/components/common/Grid/common/utils";
import { API_URL } from "src/configs/configs";
import { VIEW_MODE } from "src/common/hooks/useViewMode";

/**
 * Never expire cached collections and datasets. TODO revisit once state management approach is confirmed (#1809).
 */
const DEFAULT_QUERY_OPTIONS = {
  staleTime: Infinity,
};

/**
 * Query key for /collections/index
 */
const QUERY_ID_COLLECTIONS = "collectionsIndex";

/**
 * Query key for /datasets/index
 */
const QUERY_ID_DATASETS = "datasetIndex";

/**
 * Model returned on fetch of collections or datasets: materialized view models (rows) as well as fetch status.
 */
export interface FetchCategoriesRows<T extends Categories> {
  isError: boolean;
  isLoading: boolean;
  isSuccess: boolean;
  rows: T[];
}

/**
 * Model of collection IDs in revision keyed by corresponding published collection ID.
 */
type CollectionRevisionIdByCollectionId = Map<
  UserCollectionResponse["revision_of"],
  UserCollectionResponse["id"]
>;

/**
 * Model of /collections/index JSON response.
 */
export interface CollectionResponse {
  id: string;
  name: string;
  published_at: number;
  publisher_metadata: PublisherMetadata;
  revised_at: number;
}

/**
 * Model of /user-collections/index JSON response.
 */
export interface UserCollectionResponse extends CollectionResponse {
  access_type: ACCESS_TYPE;
  curator_name: string;
  owner: Collection["owner"];
  revision_of: Collection["revision_of"] | null;
  visibility: VISIBILITY_TYPE;
}

/**
 * Model of /collections/index or /user-collections/index JSON response that has been modified to include calculated
 * fields that facilitate filter functionality.
 */
export type ProcessedCollectionResponse = (
  | CollectionResponse
  | UserCollectionResponse
) & {
  publicationAuthors: string[];
  publicationDateValues: number[];
  revisedBy?: string;
  status?: COLLECTION_STATUS[];
};

/**
 * Model of /datasets/index JSON response.
 */
export interface DatasetResponse {
  assay: Ontology[];
  cell_count: number | null;
  cell_type: Ontology[];
  cell_type_ancestors: string[];
  collection_id: string;
  development_stage_ancestors: string[];
  disease: Ontology[];
  self_reported_ethnicity: Ontology[];
  explorer_url: string;
  id: string;
  is_primary_data: IS_PRIMARY_DATA;
  mean_genes_per_cell: number | null;
  name: string;
  organism: Ontology[];
  published_at: number;
  revised_at?: number;
  sex: Ontology[];
  suspension_type: string[];
  tissue: Ontology[];
  tissue_ancestors: string[];
}

/**
 * Model of /datasets/index JSON response that has been modified to include calculated fields that facilitate filter
 * functionality.
 */
export interface ProcessedDatasetResponse extends DatasetResponse {
  cellTypeCalculated: string[]; // Field to drive cell class, cell subclass and cell type filter functionality.
  tissueCalculated: string[]; // Field to drive tissue system, tissue organ and tissue filter functionality.
}

/**
 * Query key for caching collections returned from /collections/index endpoint.
 */
export const USE_COLLECTIONS_INDEX = {
  entities: [ENTITIES.COLLECTION],
  id: QUERY_ID_COLLECTIONS,
};

/**
 * Query key for caching datasets returned from /datasets/index endpoint.
 */
export const USE_DATASETS_INDEX = {
  entities: [ENTITIES.DATASET],
  id: QUERY_ID_DATASETS,
};

/**
 * Fetch collection and dataset information and build collection-specific filter view model.
 * @param mode - View mode.
 * @param status - Query status.
 * @returns All public collections and the aggregated metadata of their datasets.
 */
export function useFetchCollectionRows(
  mode: VIEW_MODE,
  status: QueryStatus
): FetchCategoriesRows<CollectionRow> {
  // Fetch datasets.
  const {
    data: datasets,
    isError: datasetsError,
    isLoading: datasetsLoading,
    isSuccess: datasetsSuccess,
  } = useFetchDatasets(mode, status);

  // Fetch collections.
  const {
    data: collectionsById,
    isError: collectionsError,
    isLoading: collectionsLoading,
    isSuccess: collectionsSuccess,
  } = useFetchCollections(mode, status);

  // View model built from join of collections response and aggregated metadata of dataset rows.
  // Build dataset rows once datasets and collections responses have resolved.
  const collectionRows = useMemo(() => {
    if (!datasets || !collectionsById) {
      return [];
    }
    const datasetRows = buildDatasetRows(collectionsById, datasets);
    return buildCollectionRows(collectionsById, datasetRows);
  }, [datasets, collectionsById]);

  return {
    isError: datasetsError || collectionsError,
    isLoading: datasetsLoading || collectionsLoading,
    isSuccess: datasetsSuccess && collectionsSuccess,
    rows: collectionRows,
  };
}

/**
 * Cache-enabled hook for fetching collections and returning only core collection fields.
 * @param mode - View mode.
 * @param status - Query status.
 * @returns Array of collections - possible cached from previous request - containing only ID, name and recency values.
 */
export function useFetchCollections(
  mode: VIEW_MODE = VIEW_MODE.DEFAULT,
  status: QueryStatus
): UseQueryResult<Map<string, ProcessedCollectionResponse>> {
  return useQuery<Map<string, ProcessedCollectionResponse>>(
    [USE_COLLECTIONS_INDEX],
    mode === VIEW_MODE.DEFAULT ? fetchCollections : fetchUserCollections,
    {
      ...DEFAULT_QUERY_OPTIONS,
      enabled: status === "success" || status === "error",
    }
  );
}

/**
 * Fetch collection and dataset information and build filter view model.
 * @param mode - View mode.
 * @param status - Query status.
 * @returns All public datasets joined with their corresponding collection information.
 */
export function useFetchDatasetRows(
  mode: VIEW_MODE,
  status: QueryStatus
): FetchCategoriesRows<DatasetRow> {
  // Fetch datasets.
  const {
    data: datasets,
    isError: datasetsError,
    isLoading: datasetsLoading,
    isSuccess: datasetsSuccess,
  } = useFetchDatasets(mode, status);

  // Fetch collections.
  const {
    data: collectionsById,
    isError: collectionsError,
    isLoading: collectionsLoading,
    isSuccess: collectionsSuccess,
  } = useFetchCollections(mode, status);

  // Build dataset rows once datasets and collections responses have resolved.
  const datasetRows = useMemo(() => {
    if (!datasets || !collectionsById) {
      return [];
    }
    // Default view mode displays all datasets, while the curator view mode only displays published datasets.
    // If the curator view mode is active, any private collections are removed from the collectionsById map, and
    // the corresponding datasets are not shown in the dataset rows.
    return buildDatasetRows(
      deletePrivateCollectionsById(collectionsById),
      datasets
    );
  }, [datasets, collectionsById]);

  return {
    isError: datasetsError || collectionsError,
    isLoading: datasetsLoading || collectionsLoading,
    isSuccess: datasetsSuccess && collectionsSuccess,
    rows: datasetRows,
  };
}

/**
 * Cache-enabled hook for fetching public, non-tombstoned, datasets returning only filterable and sortable fields.
 * @param mode - View mode.
 * @param status - Query status.
 * @returns Array of datasets - possible cached from previous request - containing filterable and sortable dataset
 * fields.
 */
export function useFetchDatasets(
  mode: VIEW_MODE = VIEW_MODE.DEFAULT,
  status: QueryStatus
): UseQueryResult<ProcessedDatasetResponse[]> {
  return useQuery<ProcessedDatasetResponse[]>(
    [USE_DATASETS_INDEX],
    mode === VIEW_MODE.DEFAULT ? fetchDatasets : fetchUserDatasets,
    {
      ...DEFAULT_QUERY_OPTIONS,
      enabled: status === "success" || status === "error",
    }
  );
}

/**
 * Create model of collection category values by aggregating the values in each category of each dataset in collection.
 * @param collectionDatasetRows - Datasets in the collection to aggregate category values over.
 * @returns Object containing aggregated category values from given dataset rows.
 */
function aggregateCollectionDatasetRows(
  collectionDatasetRows: DatasetRow[]
): Categories {
  // Aggregate dataset category values for each category in the collection.
  const aggregatedCategoryValues = collectionDatasetRows.reduce(
    (accum: Categories, collectionDatasetRow: DatasetRow) => {
      return {
        assay: [...accum.assay, ...collectionDatasetRow.assay],
        cell_count: accum.cell_count + collectionDatasetRow.cell_count,
        cellTypeCalculated: [
          ...accum.cellTypeCalculated,
          ...collectionDatasetRow.cellTypeCalculated,
        ],
        cell_type: [...accum.cell_type, ...collectionDatasetRow.cell_type],
        cell_type_ancestors: [
          ...accum.cell_type_ancestors,
          ...collectionDatasetRow.cell_type_ancestors,
        ],
        development_stage_ancestors: [
          ...accum.development_stage_ancestors,
          ...collectionDatasetRow.development_stage_ancestors,
        ],
        disease: [...accum.disease, ...collectionDatasetRow.disease],
        organism: [...accum.organism, ...collectionDatasetRow.organism],
        self_reported_ethnicity: [
          ...accum.self_reported_ethnicity,
          ...collectionDatasetRow.self_reported_ethnicity,
        ],
        sex: [...accum.sex, ...collectionDatasetRow.sex],
        suspension_type: [
          ...accum.suspension_type,
          ...collectionDatasetRow.suspension_type,
        ],
        tissue: [...accum.tissue, ...collectionDatasetRow.tissue],
        tissueCalculated: [
          ...accum.tissueCalculated,
          ...collectionDatasetRow.tissueCalculated,
        ],
        tissue_ancestors: [
          ...accum.tissue_ancestors,
          ...collectionDatasetRow.tissue_ancestors,
        ],
      };
    },
    {
      assay: [],
      cell_count: 0,
      cellTypeCalculated: [],
      cell_type: [],
      cell_type_ancestors: [],
      development_stage_ancestors: [],
      disease: [],
      organism: [],
      self_reported_ethnicity: [],
      sex: [],
      suspension_type: [],
      tissue: [],
      tissueCalculated: [],
      tissue_ancestors: [],
    }
  );

  // De-dupe aggregated category values.
  return {
    assay: uniqueOntologies(aggregatedCategoryValues.assay),
    cell_count: aggregatedCategoryValues.cell_count,
    cellTypeCalculated: [
      ...new Set(aggregatedCategoryValues.cellTypeCalculated),
    ],
    cell_type: uniqueOntologies(aggregatedCategoryValues.cell_type),
    cell_type_ancestors: [
      ...new Set(aggregatedCategoryValues.cell_type_ancestors),
    ],
    development_stage_ancestors: [
      ...new Set(aggregatedCategoryValues.development_stage_ancestors),
    ],
    disease: uniqueOntologies(aggregatedCategoryValues.disease),
    organism: uniqueOntologies(aggregatedCategoryValues.organism),
    self_reported_ethnicity: uniqueOntologies(
      aggregatedCategoryValues.self_reported_ethnicity
    ),
    sex: uniqueOntologies(aggregatedCategoryValues.sex),
    suspension_type: [...new Set(aggregatedCategoryValues.suspension_type)],
    tissue: uniqueOntologies(aggregatedCategoryValues.tissue),
    tissueCalculated: [...new Set(aggregatedCategoryValues.tissueCalculated)],
    tissue_ancestors: [...new Set(aggregatedCategoryValues.tissue_ancestors)],
  };
}

/**
 * Create collection rows from aggregated dataset category values and add to each dataset in collection.
 * @param collectionsById - Collections keyed by their ID.
 * @param datasetRows - Array of joined of dataset and basic collection information (that is, collection name).
 * @returns Datasets joined with their corresponding collection information as well as aggregated category values
 * across sibling datasets in its collection.
 */
function buildCollectionRows(
  collectionsById: Map<string, ProcessedCollectionResponse>,
  datasetRows: DatasetRow[]
): CollectionRow[] {
  // Group datasets by collection to facilitate aggregation of dataset category values for each collection.
  const datasetRowsByCollectionId = groupDatasetRowsByCollection(datasetRows);

  // Aggregate category values for each collection and update on each dataset.
  const collectionRows = [];
  for (const [collectionId, collection] of collectionsById.entries()) {
    // Create model of collection category values by aggregating the values in each category of each dataset in
    // collection.
    const collectionDatasetRows =
      datasetRowsByCollectionId.get(collectionId) ?? [];
    const aggregatedCategoryValues = aggregateCollectionDatasetRows(
      collectionDatasetRows
    );

    // Calculate recency for sorting.
    const recency = calculateRecency(
      collection,
      collection?.publisher_metadata
    );

    // Build the summary citation from the collection's publication metadata, if any.
    const summaryCitation = buildSummaryCitation(
      collection?.publisher_metadata
    );

    // Calculate test ID
    const testId = createCollectionRowTestId(collection);

    // Create collection row from aggregated collection category values and core collection information.
    const collectionRow = sortCategoryValues({
      ...collection,
      ...aggregatedCategoryValues,
      recency,
      summaryCitation,
      testId,
    });
    collectionRows.push(collectionRow);
  }
  return collectionRows;
}

/**
 * Join dataset and collection information to facilitate filter over datasets.
 * @param collectionsById - Collections keyed by their ID.
 * @param datasets - Datasets returned from datasets/index endpoint.
 * @returns Datasets joined with their corresponding collection information.
 */
function buildDatasetRows(
  collectionsById: Map<string, ProcessedCollectionResponse>,
  datasets: ProcessedDatasetResponse[]
): DatasetRow[] {
  // Join collection and dataset information to create dataset rows.
  return datasets.reduce((acc, dataset: ProcessedDatasetResponse) => {
    const collection = collectionsById.get(dataset.collection_id);
    if (collection) {
      acc.push(buildDatasetRow(dataset, collection));
    }
    return acc;
  }, [] as DatasetRow[]);
}

/**
 * Build dataset row from dataset response.
 * @param dataset - Response dataset values to build filterable data from.
 * @param collection - Response collection values to join with dataset values, possibly undefined if dataset is an
 * orphan with no corresponding collection.
 * @returns Fully built dataset row; join between dataset and collection values with corrected missing and
 * data primary values.
 */
function buildDatasetRow(
  dataset: ProcessedDatasetResponse,
  collection?: ProcessedCollectionResponse
): DatasetRow {
  // Determine dataset's publication month and year.
  const [publicationMonth, publicationYear] = getPublicationMonthYear(
    dataset,
    collection?.publisher_metadata
  );

  // Calculate date bins for dataset.
  const [todayMonth, todayYear] = getMonthYear(new Date());
  const publicationDateValues = expandPublicationDateValues(
    todayMonth,
    todayYear,
    publicationMonth,
    publicationYear
  );

  // Calculate recency for sorting.
  const recency = calculateRecency(dataset, collection?.publisher_metadata);

  // Join!
  const datasetRow = {
    ...dataset,
    cell_count: dataset.cell_count ?? 0,
    collection_name: collection?.name ?? "-",
    isOverMaxCellCount: checkIsOverMaxCellCount(dataset.cell_count),
    publicationAuthors: collection?.publicationAuthors,
    publicationDateValues,
    recency,
  };
  return sortCategoryValues(datasetRow);
}

/**
 * Build summary citation format from given publisher metadata:
 * Last name of first author (publication year) journal abbreviation such as Ren et al. (2021) Cell.
 * @param publisherMetadata - Publication metadata of a collection.
 */
export function buildSummaryCitation(
  publisherMetadata?: PublisherMetadata
): string {
  if (!publisherMetadata) {
    return "";
  }

  const citationTokens = [];

  // Add author to citation - family name if first author is a person, name if first author is a consortium.
  const { authors, journal, published_year: publishedYear } = publisherMetadata;
  const [firstAuthor] = authors ?? [];
  if (firstAuthor) {
    if (isAuthorPerson(firstAuthor)) {
      citationTokens.push(firstAuthor.family);
    } else {
      citationTokens.push(firstAuthor.name);
    }

    if (authors.length > 1) {
      citationTokens.push("et al.");
    }
  }

  // Add year and journal.
  citationTokens.push(`(${publishedYear})`);
  citationTokens.push(journal);

  return citationTokens.join(" ");
}

/**
 * Determine date bins that publication date falls within. Date bins are recency-based:
 * - Published 1 month ago, all date bins are applicable (1 month to 3 years).
 * - Published 2 months ago, all date bins except 1 month are applicable.
 * - Published 7 months ago, all date bins except 1, 3 and 6 months are applicable.
 */
export function createPublicationDateValues(
  monthsSincePublication: number
): number[] {
  return PUBLICATION_DATE_VALUES.reduce((accum: number[], dateBin: number) => {
    if (monthsSincePublication <= dateBin) {
      accum.push(dateBin);
    }
    return accum;
  }, []);
}

/**
 * Calculate the number of months since the collection was published. This value is used when filtering by publication
 * date.
 * @param todayMonth - Today's month
 * @param todayYear - Today's year
 * @param publicationMonth - Month of publication (1 = January, 2 = February etc).
 * @param publicationYear - Year of publication
 * @returns The count of months since collection was published.
 */
export function calculateMonthsSincePublication(
  todayMonth: number,
  todayYear: number,
  publicationMonth: number,
  publicationYear: number
): number {
  return todayMonth - publicationMonth + 12 * (todayYear - publicationYear);
}

/**
 * Create test ID for the given collection row. Test IDs are used to identify collection rows in end2end tests,
 * specially tests that require finding rows that are available for revision.
 * @param collection - Collection response returned from endpoint that has been processed to include additional
 * FE-specific values.
 * @return Test ID containing access type and status indicators for the given collection row.
 */
function createCollectionRowTestId(
  collection: ProcessedCollectionResponse
): string {
  const accessType =
    "access_type" in collection ? collection.access_type : ACCESS_TYPE.READ;
  const status = collection.status?.join("-") ?? COLLECTION_STATUS.PUBLISHED;
  return `collection-row-${accessType}-${status}`.toLowerCase();
}

/**
 * Calculate the recency value for the given collection or dataset. If there is no associated publication
 * metadata, use revised_at or published_at in priority order for the given collection or dataset.
 * @param response - Collection or dataset response returned from endpoint.
 * @param publisherMetadata - Publication metadata of the collection, or of the dataset's corresponding collection.
 * @returns Number representing seconds since Unix epoch.
 */
export function calculateRecency(
  response: CollectionResponse | DatasetResponse,
  publisherMetadata?: PublisherMetadata
): number {
  // Pull date value from publication metadata if specified.
  if (publisherMetadata) {
    return publisherMetadata.published_at;
  }

  // Collection (or dataset's collection) has no publication metadata, use revised at or published at, in priority order.
  return response.revised_at ?? response.published_at;
}

/**
 * Removes private collections from the given processed collection response.
 * @param collectionsById - Processed collection response.
 * @returns processed collection response with private collections removed.
 */
export function deletePrivateCollectionsById(
  collectionsById: Map<string, ProcessedCollectionResponse>
): Map<string, ProcessedCollectionResponse> {
  const publishedCollectionsById = new Map(collectionsById);
  for (const [id, collection] of publishedCollectionsById.entries()) {
    if (
      "visibility" in collection &&
      collection?.visibility === VISIBILITY_TYPE.PRIVATE
    ) {
      publishedCollectionsById.delete(id);
    }
  }
  return publishedCollectionsById;
}

/**
 * Concat author first and last names to facilitate filter. Ignore authors with a `name` attribute as this indicates
 * author is a consortium which are not to be included in the filter.
 * @param authors - Array of collection publication authors.
 * @returns Array of strings containing author first and last names.
 */
function expandPublicationAuthors(authors: (Author | Consortium)[]): string[] {
  return authors
    .filter(isAuthorPerson)
    .map((author: Author) => `${author.family}, ${author.given}`);
}

/**
 * Determine date bins applicable to the given publication date.
 * @param todayMonth - Today's month
 * @param todayYear - Today's year
 * @param publicationMonth - Month of publication (1 = January, 2 = February etc).
 * @param publicationYear - Year of publication
 * @returns Array of recency date bins that publication date fall within.
 */
function expandPublicationDateValues(
  todayMonth: number,
  todayYear: number,
  publicationMonth: number,
  publicationYear: number
): number[] {
  // Calculate the number of months since the publication date.
  const monthsSincePublication = calculateMonthsSincePublication(
    todayMonth,
    todayYear,
    publicationMonth,
    publicationYear
  );

  // Build array containing all date ranges that the months since publication fall within.
  return createPublicationDateValues(monthsSincePublication);
}

/**
 * Fetch and process collection response from given collection endpoint.
 * @param url - URL of collection endpoint.
 * @returns Promise that resolves to a processed collection response - possible cached from previous request.
 */
function fetchAndProcessCollectionResponse(
  url: string
): Promise<ProcessedCollectionResponse[]> {
  return fetch(url, DEFAULT_FETCH_OPTIONS)
    .then((response) => response.json())
    .then((collections: CollectionResponse[] | UserCollectionResponse[]) => {
      // Calculate the number of months since publication for each collection.
      const [todayMonth, todayYear] = getMonthYear(new Date());
      return collections.map((collection: CollectionResponse) =>
        processCollectionResponse(collection, todayMonth, todayYear)
      );
    });
}

/**
 * Fetch public collections from /datasets/index endpoint. Collections are partial in that they do not contain all
 * fields; only fields required for filtering and sorting are returned.
 * @returns Promise that resolves to a map of collections keyed by collection ID - possible cached from previous
 * request - containing only ID, name and recency values.
 */
async function fetchCollections(): Promise<
  Map<string, ProcessedCollectionResponse>
> {
  const processedCollections = await fetchAndProcessCollectionResponse(
    API_URL + API.COLLECTIONS_INDEX
  );

  // Create "collections lookup" to facilitate join between collections and datasets.
  return keyCollectionsById(processedCollections);
}

/**
 * Fetch collections from /user-collections/index endpoint. User collections are partial in that they do not contain all
 * fields; only fields required for filtering and sorting are returned.
 * @returns Promise that resolves to a map of user collections keyed by collection ID - possible cached from previous
 * request.
 */
async function fetchUserCollections(): Promise<
  Map<string, ProcessedCollectionResponse>
> {
  const processedCollections = await fetchAndProcessCollectionResponse(
    API_URL + API.USER_COLLECTIONS_INDEX
  );

  // Map collection IDs in revision keyed by corresponding published collection ID.
  const collectionRevisionIdById: CollectionRevisionIdByCollectionId =
    new Map();

  // Published collection in revision are filtered out from the processed response, and the revision ID is mapped to
  // the corresponding published ID.
  // Process user collections response to add additional fields required for filtering and sorting .e.g "status".
  const processedUserCollections = processedCollections
    .filter((processedCollection) =>
      filterCollectionInRevision(processedCollection, collectionRevisionIdById)
    )
    .map((processedCollection) =>
      processUserCollectionResponse(
        processedCollection,
        collectionRevisionIdById
      )
    );

  // Create "collections lookup" to facilitate join between user collections and datasets.
  return keyCollectionsById(processedUserCollections);
}

/**
 * Fetch and process dataset response from given dataset endpoint.
 * @param url - URL of dataset endpoint.
 * @returns Promise that resolves to a processed dataset response - possible cached from previous request.
 */
function fetchAndProcessDatasetResponse(
  url: string
): Promise<ProcessedDatasetResponse[]> {
  return fetch(url, DEFAULT_FETCH_OPTIONS)
    .then((response) => response.json())
    .then((datasets: DatasetResponse[]) => {
      // Correct any dirty data returned from endpoint.
      const sanitizedDatasets = datasets.map((dataset: DatasetResponse) => {
        return sanitizeDatasetResponse(dataset);
      });
      // Process dataset response.
      return sanitizedDatasets.map((dataset: DatasetResponse) =>
        processDatasetResponse(dataset)
      );
    });
}

/**
 * Fetch public, non-tombstoned, partial datasets from /datasets/index endpoint. Datasets are partial in that they
 * do not contain all fields; only fields required for filtering and sorting are returned. Correct any dirt data
 * returned from endpoint.
 * @returns Promise resolving to an array of datasets - possible cached from previous request - containing
 * filterable and sortable dataset fields.
 */
async function fetchDatasets(): Promise<ProcessedDatasetResponse[]> {
  return await fetchAndProcessDatasetResponse(API_URL + API.DATASETS_INDEX);
}

/**
 * Fetch public, non-tombstoned, partial datasets from /user-datasets/index endpoint. Datasets are partial in that they
 * do not contain all fields; only fields required for filtering and sorting are returned. Correct any dirt data
 * returned from endpoint.
 * @returns Promise resolving to an array of datasets - possible cached from previous request - containing
 * filterable and sortable dataset fields.
 */
async function fetchUserDatasets(): Promise<ProcessedDatasetResponse[]> {
  return await fetchAndProcessDatasetResponse(
    API_URL + API.USER_DATASETS_INDEX
  );
}

/**
 * Filters collection in revision from the processed response.
 * Maps a collection in revision, by ID, to the corresponding published collection ID.
 * @param processedCollection - Processed collection.
 * @param collectionRevisionIdById - Map of collection IDs in revision keyed by corresponding published collection ID.
 * @returns True if collection is not in revision.
 */
function filterCollectionInRevision(
  processedCollection: ProcessedCollectionResponse,
  collectionRevisionIdById: CollectionRevisionIdByCollectionId
) {
  if ("revision_of" in processedCollection) {
    const { id, revision_of } = processedCollection;
    if (revision_of && revision_of !== id) {
      collectionRevisionIdById.set(revision_of, id);
      return false; // Collection is in revision, filter out.
    }
  }
  return true;
}

/**
 * Returns the status of the given collection.
 * @param processedCollection - Processed collection.
 * @param revisedBy - ID of the corresponding collection in revision.
 * @returns Collection statuses.
 */
function getCollectionStatus(
  processedCollection: ProcessedCollectionResponse,
  revisedBy?: string
): COLLECTION_STATUS[] | undefined {
  if ("access_type" in processedCollection) {
    const status: COLLECTION_STATUS[] = [];
    if (processedCollection.access_type === ACCESS_TYPE.READ) {
      status.push(COLLECTION_STATUS.PUBLISHED);
    } else {
      if (processedCollection.visibility === VISIBILITY_TYPE.PRIVATE) {
        status.push(COLLECTION_STATUS.PRIVATE);
      } else {
        status.push(COLLECTION_STATUS.PUBLISHED);
        revisedBy && status.push(COLLECTION_STATUS.REVISION);
      }
    }
    return status;
  }
}

/**
 * Return the publication month and year for the given collection or dataset. If there is no associated publication
 * metadata, use revised_at or published_at in priority order for the given collection or dataset.
 * @param response - Collection or dataset response returned from endpoint.
 * @param publisherMetadata - Publication metadata of the collection, or of the dataset's corresponding collection.
 * @returns Tuple containing collections publication month (1-indexed) and year.
 */
function getPublicationMonthYear(
  response: CollectionResponse | DatasetResponse,
  publisherMetadata?: PublisherMetadata
): [number, number] {
  // Pull month and year from publication metadata if specified.
  if (publisherMetadata) {
    return [
      publisherMetadata.published_month,
      publisherMetadata.published_year,
    ];
  }

  // Collection (or dataset's collection) has no publication metadata, use revised at or published at, in priority order.
  const recency = new Date(
    (response.revised_at ?? response.published_at) * 1000
  );
  return getMonthYear(recency);
}

/**
 * Return the month and year of the specified date.
 * @param date - Date to return UTC month and year of.
 * @returns Tuple containing the given date's month (1-indexed) and year.
 */
function getMonthYear(date: Date): [number, number] {
  return [date.getUTCMonth() + 1, date.getUTCFullYear()]; // JS dates are 0-indexed, publication dates are 1-indexed.
}

/**
 * Group dataset rows by collection.
 * @param datasetRows - Array of dataset rows to group by their collection ID.
 * @returns Dataset rows keyed by their collection IDs.
 */
function groupDatasetRowsByCollection(
  datasetRows: DatasetRow[]
): Map<string, DatasetRow[]> {
  return datasetRows.reduce((accum: Map<string, DatasetRow[]>, datasetRow) => {
    const datasetsByCollectionId = accum.get(datasetRow.collection_id);
    if (datasetsByCollectionId) {
      datasetsByCollectionId.push(datasetRow);
    } else {
      accum.set(datasetRow.collection_id, [datasetRow]);
    }
    return accum;
  }, new Map<string, DatasetRow[]>());
}

/**
 * Publication authors can be a person or a consortium; determine if the given author is in fact a person (and not a
 * consortium).
 * @param author - Person or consortium associated with a publication.
 * @returns True if author is a person and not a consortium.
 */
function isAuthorPerson(author: Author | Consortium): author is Author {
  return (author as Author).given !== undefined;
}

/**
 * Created a map of collections keyed by their ID.
 * @param collections - Collections returned from collection/index endpoint.
 * @returns Map of collections keyed by their ID.
 */
function keyCollectionsById(
  collections: ProcessedCollectionResponse[]
): Map<string, ProcessedCollectionResponse> {
  return new Map(
    collections.map((collection: ProcessedCollectionResponse) => [
      collection.id,
      collection,
    ])
  );
}

/**
 * Add calculated fields to collection response.
 * @param collection - Collection response returned from collections or user-collections endpoint.
 * @param todayMonth - Number indicating the month of today's date (1-indexed).
 * @param todayYear - Number indicating the year of today's date.
 */
function processCollectionResponse(
  collection: CollectionResponse | UserCollectionResponse,
  todayMonth: number,
  todayYear: number
): ProcessedCollectionResponse {
  // Determine the collection's publication month and year.
  const [publicationMonth, publicationYear] = getPublicationMonthYear(
    collection,
    collection.publisher_metadata
  );

  // Calculate date bins and add to "processed" collection model.
  const publicationDateValues = expandPublicationDateValues(
    todayMonth,
    todayYear,
    publicationMonth,
    publicationYear
  );

  // Determine the set of authors of the publication.
  const publicationAuthors = expandPublicationAuthors(
    collection?.publisher_metadata?.authors ?? []
  );

  return {
    ...collection,
    publicationAuthors,
    publicationDateValues,
  };
}

/**
 * Add calculated fields to dataset response.
 * @param dataset - Dataset response returned from collections endpoint.
 */
function processDatasetResponse(
  dataset: DatasetResponse
): ProcessedDatasetResponse {
  // Build up values to facilitate ontology-aware cel type and tissue filtering.
  const cellTypeCalculated = [
    ...tagAncestorsAsInferred(dataset.cell_type_ancestors),
    ...tagOntologyTermsAsExplicit(dataset.cell_type),
  ];
  const tissueCalculated = [
    ...tagAncestorsAsInferred(dataset.tissue_ancestors),
    ...tagOntologyTermsAsExplicit(dataset.tissue),
  ];
  return {
    ...dataset,
    cellTypeCalculated,
    tissueCalculated,
  };
}

/**
 * Add calculated fields to user collections response.
 * @param processedCollection - Processed collection response.
 * @param collectionRevisionIdById - Map of collection IDs in revision keyed by corresponding published collection ID.
 */
function processUserCollectionResponse(
  processedCollection: ProcessedCollectionResponse,
  collectionRevisionIdById: CollectionRevisionIdByCollectionId
): ProcessedCollectionResponse {
  // Determine the collection's corresponding revision ID.
  const revisedBy = collectionRevisionIdById.get(processedCollection.id);

  // Determine the collection's status.
  const status = getCollectionStatus(processedCollection, revisedBy);

  return { ...processedCollection, revisedBy, status };
}

/**
 * Add defaults for missing filterable values, e.g. convert missing ontology values to empty array.
 * Remove any self-reported ethnicity values on the deny list.
 * @param datasetResponse - Dataset to check for missing values.
 * @returns Corrected dataset response.
 */
function sanitizeDatasetResponse(
  datasetResponse: DatasetResponse
): DatasetResponse {
  const sanitizedDatasetResponse = { ...datasetResponse };

  sanitizedDatasetResponse.self_reported_ethnicity = (
    datasetResponse.self_reported_ethnicity ?? []
  ).filter(
    (self_reported_ethnicity) =>
      !SELF_REPORTED_ETHNICITY_DENY_LIST.includes(self_reported_ethnicity.label)
  );

  sanitizedDatasetResponse.suspension_type = (
    datasetResponse.suspension_type ?? []
  ).filter(
    (suspensionType) => !SUSPENSION_TYPE_DENY_LIST.includes(suspensionType)
  );

  sanitizedDatasetResponse.assay = datasetResponse.assay ?? [];

  sanitizedDatasetResponse.cell_type = datasetResponse.cell_type ?? [];

  sanitizedDatasetResponse.cell_type_ancestors =
    datasetResponse.cell_type_ancestors ?? [];

  sanitizedDatasetResponse.development_stage_ancestors =
    datasetResponse.development_stage_ancestors ?? [];

  sanitizedDatasetResponse.disease = datasetResponse.disease ?? [];

  sanitizedDatasetResponse.organism = datasetResponse.organism ?? [];

  sanitizedDatasetResponse.sex = datasetResponse.sex ?? [];

  sanitizedDatasetResponse.tissue = datasetResponse.tissue ?? [];

  sanitizedDatasetResponse.tissue_ancestors =
    datasetResponse.tissue_ancestors ?? [];

  return sanitizedDatasetResponse;
}

/**
 * Sort category values on the given collection or dataset rows.
 * @param row - Collection or dataset row to sort category values of.
 * @returns Array of collection or dataset rows with category values sorted.
 */
function sortCategoryValues<T extends Categories>(row: T): T {
  return {
    ...row,
    assay: row.assay.sort(sortOntologies),
    cell_type: row.cell_type.sort(sortOntologies),
    disease: row.disease.sort(sortOntologies),
    organism: row.organism.sort(sortOntologies),
    sex: row.sex.sort(sortOntologies),
    tissue: row.tissue.sort(sortOntologies),
  };
}

/*
 * Sort ontologies by label, case insensitive, ascending.
 * @param o0 - First filtered rows to compare.
 * @param o1 - Second filtered rows to compare.
 * @returns Number indicating sort precedence of o0 vs o1.
 */
function sortOntologies(o0: Ontology, o1: Ontology): number {
  return COLLATOR_CASE_INSENSITIVE.compare(o0.label, o1.label);
}

/**
 * Convert each ontology term ID into one marked as "inferred", use by tissue multi-panel category filter. For example,
 * both the tissue system and the tissue organ panels contain inferred values and allow us to execute an "includes"
 * filter in tissue ancestors rather than an exact match. So, selecting blood (in the tissue organ panel) returns
 * values that are exactly blood, or that have blood as an ancestor.
 * @param ontologyTermIds - Array of ontology term IDs to convert to an "inferred" format.
 * @returns An array of ontology terms prefixed with the "inferred" indicator.
 */
function tagAncestorsAsInferred(ontologyTermIds: string[]): string[] {
  return ontologyTermIds.map((ontologyTermId) =>
    buildInferredOntologyTermId(ontologyTermId)
  );
}

/**
 * Convert each ontology term ID into one marked as "explicit", use by tissue multi-panel category filter. For example,
 * the tissue panel contains explicit values and allow us to execute an exact filter. So, selecting blood non specific
 * (in the tissue panel) returns only values that are exactly blood.
 * @param ontologyTermIds - Array of ontology term IDs to convert to an "explicit" format.
 * @returns An array of ontology terms prefixed with the "inferred" indicator.
 */
function tagOntologyTermsAsExplicit(ontologyTermIds: Ontology[]): string[] {
  return ontologyTermIds.map((ontologyTermId) =>
    buildExplicitOntologyTermId(ontologyTermId.ontology_term_id)
  );
}

/**
 * De-dupe ontologies in the given array.
 * @param ontologies - Array of ontologies to remove duplicated from.
 * @returns Array containing set of ontologies.
 */
function uniqueOntologies(ontologies: Ontology[]): Ontology[] {
  return [
    ...new Map(
      ontologies.map((ontology: Ontology) => [ontology.label, ontology])
    ).values(),
  ];
}
