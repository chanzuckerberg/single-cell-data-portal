import { Collection, IS_PRIMARY_DATA } from "src/common/entities";
import {
  CategoryKey,
  CATEGORY_KEY,
  CollectionRow,
  DatasetRow,
} from "src/components/common/Filter/common/entities";
import collectionsIndex from "../../../tests/features/fixtures/collections/collections-index";
import datasetsIndex from "../../../tests/features/fixtures/datasets/datasets-index";

/* Model of /collections/index JSON response  */
export type CollectionsResponse = Pick<
  Collection,
  "id" | "name" | "published_at" | "revised_at"
>;

/* Model of /datasets/index JSON response  */
export type DatasetsResponse =
  | Omit<
      DatasetRow,
      | "collection_name"
      | "collection_published_at"
      | "collection_revised_at"
      | "is_primary_data"
    > & { is_primary_data: IS_PRIMARY_DATA }; // is_primary_data returned as a single value, map to array for filtering.

/**
 * Fetch collection and dataset information and build collection-specific filter view model.
 * @returns All public datasets joined with their corresponding collection information.
 */
// TODO rename to useFetchCollectionRows() and convert to react-query (#1609)
export function fetchCollectionRows(): CollectionRow[] {
  return buildCollectionRows(collectionsIndex, datasetsIndex);
}

/**
 * Fetch collection and dataset information and build filter view model.
 * @returns All public datasets joined with their corresponding collection information.
 */
// TODO rename to useFetchDatasetRows() and convert to react-query (#1609)
export function fetchDatasetRows(): DatasetRow[] {
  return buildDatasetRows(collectionsIndex, datasetsIndex);
}

/**
 * Determine the set of category values that exist in the given datasets for the given category.
 * @param categoryKey - Key of the category to aggregate values of.
 * @param datasetRows - Datasets to aggregate category values of.
 * @returns Array of aggregated category values for the given category.
 * TODO(cc) set type and return type
 */
function aggregateCategory(
  categoryKey: CategoryKey,
  datasetRows: DatasetRow[]
) {
  const categoryValuesSet = datasetRows.reduce(
    (accum, datasetRow: DatasetRow) => {
      datasetRow[categoryKey].forEach((categoryValue) =>
        accum.add(categoryValue)
      );
      return accum;
    },
    new Set()
  );
  return [...categoryValuesSet];
}

/**
 * Create model of collection category values by aggregating the values in each category of each dataset in collection.
 * @param datasetRows - Datasets in the collection to aggregate category values over.
 * @returns Object containing aggregated category values from given dataset rows.
 */
function aggregateCollectionCategoryValues(
  datasetRows: DatasetRow[]
): Partial<CollectionRow> {
  // Aggregate dataset category values for each category in the collection.
  return Object.values(CATEGORY_KEY).reduce(
    (accum: Partial<CollectionRow>, categoryKey: CategoryKey) => {
      // @ts-expect-error -- TODO(cc) revisit
      accum[`${categoryKey}Aggregated`] = aggregateCategory(
        categoryKey,
        datasetRows
      );
      return accum;
    },
    {}
  );
}

/**
 * Create collection rows from aggregated dataset category values and add to each dataset in collection.
 * @param responseCollections - Collections returned from collection/index endpoint.
 * @param responseDatasets - Datasets returned from datasets/index endpoint.
 * @returns Datasets joined with their corresponding collection information as well as aggregated category values
 * across sibling datasets in its collection.
 */
function buildCollectionRows(
  responseCollections: CollectionsResponse[],
  responseDatasets: DatasetsResponse[]
): CollectionRow[] {
  // Build the base filtering model; aggregated dataset category values will be added to these bases.
  const baseDatasetRows = buildDatasetRows(
    responseCollections,
    responseDatasets
  );

  // Group datasets by collection to facilitate aggregation of dataset category values per collection.
  const datasetsByCollectionId = groupDatasetsByCollection(baseDatasetRows);

  // Aggregate category values for each collection and update on each dataset.
  const groupedDatasetRows = [...datasetsByCollectionId.values()].map(
    (datasetRows: DatasetRow[]) => {
      // Create model of collection category values by aggregating the values in each category of each dataset in
      // collection.
      const aggregatedCollectionCategoryValues =
        aggregateCollectionCategoryValues(datasetRows);

      // Add aggregated collection category values to each dataset
      return datasetRows.map((datasetRow: DatasetRow) => {
        return {
          ...datasetRow,
          ...aggregatedCollectionCategoryValues,
        };
      });
    }
  );

  // Flatten the array of dataset rows array.
  // @ts-expect-error -- TODO(cc) revisit
  return groupedDatasetRows.flat();
}

/**
 * Join dataset and collection information to facilitate filter over datasets.
 * @param responseCollections - Collections returned from collection/index endpoint.
 * @param responseDatasets - Datasets returned from datasets/index endpoint.
 * @returns Datasets joined with their corresponding collection information.
 */
function buildDatasetRows(
  responseCollections: CollectionsResponse[],
  responseDatasets: DatasetsResponse[]
): DatasetRow[] {
  // Group collections by ID to facilitate join with datasets.
  const collectionsById = groupCollectionsById(responseCollections);

  // Join collection and dataset information to create dataset row.
  return responseDatasets.reduce(
    (accum: DatasetRow[], responseDataset: DatasetsResponse) => {
      // Grab the collection for this dataset.
      const responseCollection = collectionsById.get(
        responseDataset.collection_id
      );
      if (!responseCollection) {
        console.error(
          `Collection "${responseDataset.collection_id}" not found for dataset "${responseDataset.id}"`
        );
        return accum;
      }

      accum.push(buildFilterableDataset(responseDataset, responseCollection));

      return accum;
    },
    []
  );
}

/**
 * Build dataset row from datasets index response. Add missing values and correct primary data flag where
 * necessary.
 * @param responseDataset - Response dataset values to build filterable data from.
 * @param responseCollection - Response collection values to join with dataset values.
 * @returns Fully built dataset row; join between dataset and collection values with corrected missing and
 * data primary values.
 */
function buildFilterableDataset(
  responseDataset: DatasetsResponse,
  responseCollection: CollectionsResponse
): DatasetRow {
  // Add any missing values.
  const { is_primary_data, ...sanitizedResponseDataset } =
    sanitizeDatasetIndexResponse(responseDataset);

  // Correct is_primary_data values.
  const isPrimaryData = sanitizeIsPrimaryData(is_primary_data);

  // Join!
  return {
    ...sanitizedResponseDataset,
    collection_name: responseCollection.name,
    collection_published_at: responseCollection.published_at,
    collection_revised_at: responseCollection.revised_at,
    is_primary_data: isPrimaryData,
  };
}

/**
 * Group collection responses by ID.
 * @param responseCollections - Collections returned from collection/index endpoint.
 */
function groupCollectionsById(
  responseCollections: CollectionsResponse[]
): Map<string, CollectionsResponse> {
  return responseCollections.reduce(
    (
      accum: Map<string, CollectionsResponse>,
      responseCollection: CollectionsResponse
    ) => accum.set(responseCollection.id, responseCollection),
    new Map<string, CollectionsResponse>()
  );
}

/**
 * Group dataset row by collection.
 * @param datasetRows - Array of dataset rows to group by their collection ID.
 * @returns Map of dataset rows key by collection ID.
 */
function groupDatasetsByCollection(
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
 * Add defaults for any missing filterable values on the given dataset.
 * @param responseDataset - Dataset to check for missing values.
 * @returns Corrected dataset response.
 */
function sanitizeDatasetIndexResponse(
  responseDataset: DatasetsResponse
): DatasetsResponse {
  return Object.values(CATEGORY_KEY).reduce(
    (accum: DatasetsResponse, categoryKey: CATEGORY_KEY) => {
      if (categoryKey === CATEGORY_KEY.IS_PRIMARY_DATA) {
        return accum;
      }
      accum[categoryKey] = responseDataset[categoryKey] ?? [];
      return accum;
    },
    { ...responseDataset }
  );
}

/**
 * Determine the correct value for is_primary_data. Convert "BOTH" primary data values to ["primary", "secondary"].
 * Convert "primary" or "secondary" to singleton arrays. Convert error cases where is_primary_data is undefined to [].
 * @param isPrimaryData - Primary data value to sanitize.
 */
function sanitizeIsPrimaryData(
  isPrimaryData: IS_PRIMARY_DATA
): IS_PRIMARY_DATA[] {
  if (!isPrimaryData) {
    return [];
  }

  return isPrimaryData === IS_PRIMARY_DATA.BOTH
    ? [IS_PRIMARY_DATA.PRIMARY, IS_PRIMARY_DATA.SECONDARY]
    : [isPrimaryData];
}
