import { IS_PRIMARY_DATA, Ontology } from "src/common/entities";
import {
  CATEGORY_KEY,
  CollectionCategoryValues,
  CollectionRow,
  DatasetRow,
} from "src/components/common/Filter/common/entities";
import collectionsIndex from "../../../tests/features/fixtures/collections/collections-index";
import datasetsIndex from "../../../tests/features/fixtures/datasets/datasets-index";

/* Model of /collections/index JSON response  */
export interface CollectionResponse {
  id: string;
  name: string;
  published_at: number;
  revised_at: number;
}

/* Model of /datasets/index JSON response  */
export interface DatasetResponse {
  assay: Ontology[];
  cell_count: number | null;
  cell_type: Ontology[];
  collection_id: string;
  disease: Ontology[];
  id: string;
  is_primary_data: IS_PRIMARY_DATA;
  name: string;
  organism: Ontology[];
  published_at: number;
  revised_at?: number;
  sex: Ontology[];
  tissue: Ontology[];
}

/**
 * Fetch collection and dataset information and build collection-specific filter view model.
 * @returns All public datasets joined with their corresponding collection information.
 */
export function fetchCollectionRows(): CollectionRow[] {
  const collectionsById = fetchCollections();
  const datasets = fetchDatasets();
  const datasetsRows = buildDatasetRows(collectionsById, datasets);
  return buildCollectionRows(collectionsById, datasetsRows);
}

/**
 * Fetch collection and dataset information and build filter view model.
 * @returns All public datasets joined with their corresponding collection information.
 */
export function fetchDatasetRows(): DatasetRow[] {
  const collectionsById = fetchCollections();
  const datasets = fetchDatasets();
  return buildDatasetRows(collectionsById, datasets);
}

/**
 * Fetch datasets from /datasets/index. Correct any dirt data returned from endpoint.
 * @returns Array of datasets.
 */
function fetchDatasets(): DatasetResponse[] {
  // Correct any dirty data returned from endpoint.
  return datasetsIndex.map((dataset: DatasetResponse) => {
    return sanitizeDataset(dataset);
  });
}

/**
 * Fetch collections from /datasets/index. Key collections by their ID for easy lookups during join and aggregate
 * functions.
 * @returns Map of collections keyed by their ID.
 */
function fetchCollections(): Map<string, CollectionResponse> {
  // Create "collections lookup" to facilitate join between collections and datasets.
  return keyCollectionsById(collectionsIndex);
}

/**
 * Create model of collection category values by aggregating the values in each category of each dataset in collection.
 * @param collectionDatasetRows - Datasets in the collection to aggregate category values over.
 * @returns Object containing aggregated category values from given dataset rows.
 */
function aggregateCollectionDatasetRows(
  collectionDatasetRows: DatasetRow[]
): CollectionCategoryValues {
  // Aggregate dataset category values for each category in the collection.
  return collectionDatasetRows.reduce(
    (accum: CollectionCategoryValues, collectionDatasetRow: DatasetRow) => {
      return {
        assayAggregated: [
          ...accum.assayAggregated,
          ...collectionDatasetRow.assay,
        ],
        cellTypeAggregated: [
          ...accum.cellTypeAggregated,
          ...collectionDatasetRow.cell_type,
        ],
        diseaseAggregated: [
          ...accum.diseaseAggregated,
          ...collectionDatasetRow.disease,
        ],
        isPrimaryDataAggregated: [
          ...accum.isPrimaryDataAggregated,
          ...collectionDatasetRow.is_primary_data,
        ],
        organismAggregated: [
          ...accum.organismAggregated,
          ...collectionDatasetRow.organism,
        ],
        sexAggregated: [...accum.sexAggregated, ...collectionDatasetRow.sex],
        tissueAggregated: [
          ...accum.tissueAggregated,
          ...collectionDatasetRow.tissue,
        ],
      };
    },
    {
      assayAggregated: [],
      cellTypeAggregated: [],
      diseaseAggregated: [],
      isPrimaryDataAggregated: [],
      organismAggregated: [],
      sexAggregated: [],
      tissueAggregated: [],
    }
  );
}

/**
 * Create collection rows from aggregated dataset category values and add to each dataset in collection.
 * @param collectionsById - Collections keyed by their ID.
 * @param datasetRows - Array of joined of dataset and basic collection information (that is, collection name).
 * @returns Datasets joined with their corresponding collection information as well as aggregated category values
 * across sibling datasets in its collection.
 */
function buildCollectionRows(
  collectionsById: Map<string, CollectionResponse>,
  datasetRows: DatasetRow[]
): CollectionRow[] {
  // Group datasets by collection to facilitate aggregation of dataset category values for each collection.
  const collectionDatasetRows = groupDatasetRowsByCollection(datasetRows);

  // Aggregate category values for each collection and update on each dataset.
  return collectionDatasetRows
    .map((datasetRows: DatasetRow[]) => {
      // Create model of collection category values by aggregating the values in each category of each dataset in
      // collection.
      const aggregatedCategoryValues =
        aggregateCollectionDatasetRows(datasetRows);

      // Add aggregated collection category values to each dataset
      return datasetRows.map((datasetRow: DatasetRow) => {
        const collection = collectionsById.get(datasetRow.collection_id);
        return {
          ...datasetRow,
          published_at: collection?.published_at,
          revised_at: collection?.revised_at,
          ...aggregatedCategoryValues,
        };
      });
    })
    .flat() as CollectionRow[];
}

/**
 * Join dataset and collection information to facilitate filter over datasets.
 * @param collectionsById - Collections keyed by their ID.
 * @param datasets - Datasets returned from datasets/index endpoint.
 * @returns Datasets joined with their corresponding collection information.
 */
function buildDatasetRows(
  collectionsById: Map<string, CollectionResponse>,
  datasets: DatasetResponse[]
): DatasetRow[] {
  // Join collection and dataset information to create dataset rows.
  return datasets.map((dataset: DatasetResponse) => {
    const collection = collectionsById.get(dataset.collection_id);
    return buildDatasetRow(dataset, collection);
  });
}

/**
 * Build dataset row from dataset response. Correct is_primary_data where necessary.
 * @param dataset - Response dataset values to build filterable data from.
 * @param collection - Response collection values to join with dataset values, possibly undefined if dataset is an
 * orphan with no corresponding collection.
 * @returns Fully built dataset row; join between dataset and collection values with corrected missing and
 * data primary values.
 */
function buildDatasetRow(
  dataset: DatasetResponse,
  collection?: CollectionResponse
): DatasetRow {
  const { is_primary_data } = dataset;

  // Join!
  return {
    ...dataset,
    collection_name: collection?.name ?? "-",
    is_primary_data: expandIsPrimaryData(is_primary_data),
  };
}

/**
 * Created a map of collections keyed by their ID.
 * @param responseCollections - Collections returned from collection/index endpoint.
 * @returns Map of collections keyed by their ID.
 */
function keyCollectionsById(
  responseCollections: CollectionResponse[]
): Map<string, CollectionResponse> {
  return responseCollections.reduce(
    (
      accum: Map<string, CollectionResponse>,
      responseCollection: CollectionResponse
    ) => accum.set(responseCollection.id, responseCollection),
    new Map<string, CollectionResponse>()
  );
}

/**
 * Group dataset rows by collection.
 * @param datasetRows - Array of dataset rows to group by their collection ID.
 * @returns Array of dataset row arrays, one array for each collection.
 */
function groupDatasetRowsByCollection(
  datasetRows: DatasetRow[]
): DatasetRow[][] {
  const datasetRowsByCollectionId = datasetRows.reduce(
    (accum: Map<string, DatasetRow[]>, datasetRow) => {
      const datasetsByCollectionId = accum.get(datasetRow.collection_id);
      if (datasetsByCollectionId) {
        datasetsByCollectionId.push(datasetRow);
      } else {
        accum.set(datasetRow.collection_id, [datasetRow]);
      }
      return accum;
    },
    new Map<string, DatasetRow[]>()
  );

  return [...datasetRowsByCollectionId.values()];
}

/**
 * Add defaults for missing filterable values: convert missing ontology values to empty array and is_primary_data to "".
 * @param dataset - Dataset to check for missing values.
 * @returns Corrected dataset response.
 */
function sanitizeDataset(dataset: DatasetResponse): DatasetResponse {
  return Object.values(CATEGORY_KEY).reduce(
    (accum: DatasetResponse, categoryKey: CATEGORY_KEY) => {
      if (categoryKey === CATEGORY_KEY.IS_PRIMARY_DATA) {
        accum.is_primary_data = dataset.is_primary_data ?? "";
      } else {
        accum[categoryKey] = dataset[categoryKey] ?? [];
      }
      return accum;
    },
    { ...dataset }
  );
}

/**
 * Determine the correct value for is_primary_data. Convert "BOTH" primary data values to ["primary", "secondary"].
 * Convert "primary" or "secondary" to singleton arrays. Convert error cases where is_primary_data is undefined to [].
 * @param isPrimaryData - Primary data value to sanitize.
 */
function expandIsPrimaryData(
  isPrimaryData: IS_PRIMARY_DATA
): IS_PRIMARY_DATA[] {
  if (!isPrimaryData) {
    return [];
  }

  return isPrimaryData === IS_PRIMARY_DATA.BOTH
    ? [IS_PRIMARY_DATA.PRIMARY, IS_PRIMARY_DATA.SECONDARY]
    : [isPrimaryData];
}
