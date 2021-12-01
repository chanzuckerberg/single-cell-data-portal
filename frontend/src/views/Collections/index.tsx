import Head from "next/head";
import React, { useMemo } from "react";
import {
  AggregatorFn,
  CellProps,
  CellValue,
  Column,
  Row,
  useFilters,
  useGroupBy,
  useSortBy,
  useTable,
} from "react-table";
import { ROUTES } from "src/common/constants/routes";
import { IS_PRIMARY_DATA, Ontology } from "src/common/entities";
import { FEATURES } from "src/common/featureFlags/features";
import {
  CategoryKey,
  CATEGORY_KEY,
  NonOntologyCategoryKey,
  OntologyCategoryKey,
  useFacetedFilter,
} from "src/common/hooks/useFacetedFilter";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import Categories from "src/components/Categories";
import FilteredCollectionsGrid from "src/components/Collections/components/Grid/components/FilteredCollectionsGrid";
import {
  FilterableCollection,
  FilterableDataset,
} from "src/components/common/Filter/common/entities";
import SideBar from "src/components/common/SideBar";
import { View } from "src/views/globalStyle";
import filterableDatasets from "../../../tests/features/fixtures/datasets/filterable-datasets";

// Collection ID object key
const COLLECTION_ID = "collection_id";

// Collection name object key
const COLLECTION_NAME = "collection_name";

// Key identifying recency sort by column
const COLUMN_ID_RECENCY = "recency";

export default function Collections(): JSX.Element {
  // Column configuration backing table.
  const columnConfig: Column<FilterableDataset>[] = useMemo(
    () => [
      // Hidden, required for grouping datasets by collections
      {
        accessor: COLLECTION_ID,
      },
      // Hidden, required for sorting TODO(cc) this should be for materialized collection row
      {
        accessor: (dataset: FilterableDataset): number =>
          dataset.revised_at ?? dataset.published_at,
        id: COLUMN_ID_RECENCY,
      },
      // Collection name, aggregated across datasets to roll up into single value for each collection ID group.
      {
        Header: "Collection",
        accessor: COLLECTION_NAME,
        aggregate: aggregateFn(),
      },
      {
        Cell: Cell,
        Header: "Tissue",
        accessor: aggregatedOntologyCellAccessorFn(CATEGORY_KEY.TISSUE),
        aggregate: aggregateFn(),
        filter: "includesSome",
        id: `${CATEGORY_KEY.TISSUE}`,
      },
      {
        Cell: Cell,
        Header: "Disease",
        accessor: aggregatedOntologyCellAccessorFn(CATEGORY_KEY.DISEASE),
        aggregate: aggregateFn(),
        filter: "includesSome",
        id: `${CATEGORY_KEY.DISEASE}`,
      },
      {
        Cell: Cell,
        Header: "Assay",
        accessor: aggregatedOntologyCellAccessorFn(CATEGORY_KEY.ASSAY),
        aggregate: aggregateFn(),
        filter: "includesSome",
        id: `${CATEGORY_KEY.ASSAY}`,
      },
      {
        Cell: Cell,
        Header: "Organism",
        accessor: aggregatedOntologyCellAccessorFn(CATEGORY_KEY.ORGANISM),
        aggregate: aggregateFn(),
        filter: "includesSome",
        id: `${CATEGORY_KEY.ORGANISM}`,
      },
      {
        Cell: Cell,
        Header: "Cell Type",
        accessor: aggregatedOntologyCellAccessorFn(CATEGORY_KEY.CELL_TYPE),
        aggregate: aggregateFn(),
        filter: "includesSome",
        id: `${CATEGORY_KEY.CELL_TYPE}`,
      },
      {
        Cell: Cell,
        Header: "Primary Data",
        accessor: aggregatedCellAccessorFn(CATEGORY_KEY.IS_PRIMARY_DATA),
        aggregate: aggregateFn(),
        filter: "includesSome",
        id: `${CATEGORY_KEY.IS_PRIMARY_DATA}`,
      },
      {
        Cell: Cell,
        Header: "Sex",
        accessor: aggregatedOntologyCellAccessorFn(CATEGORY_KEY.SEX),
        aggregate: aggregateFn(),
        filter: "includesSome",
        id: `${CATEGORY_KEY.SEX}`,
      },
    ],
    []
  );

  // Init collection-based filterable datasets.
  const data: FilterableDataset[] = useMemo(
    () => prepareData(filterableDatasets),
    []
  );

  // Table init
  const tableInstance = useTable<FilterableDataset>(
    {
      columns: columnConfig,
      data,
      initialState: {
        groupBy: [COLLECTION_ID],
        // Only display aggregated tissue, disease and organism values.
        hiddenColumns: [
          COLLECTION_ID,
          COLUMN_ID_RECENCY,
          // CATEGORY_KEY.ASSAY,
          // CATEGORY_KEY.CELL_TYPE,
          // CATEGORY_KEY.IS_PRIMARY_DATA,
          // CATEGORY_KEY.SEX,
        ],
        sortBy: [
          {
            desc: true,
            id: COLUMN_ID_RECENCY,
          },
        ],
      },
    },
    useFilters,
    useGroupBy,
    useSortBy
  );

  // Filter init.
  const {
    preFilteredRows,
    setFilter,
    state: { filters },
  } = tableInstance;
  const filterInstance = useFacetedFilter(
    preFilteredRows,
    filters,
    setFilter,
    COLLECTION_ID
  );

  // Hide datasets behind feature flag - start
  const isFilterEnabled = useFeatureFlag(FEATURES.FILTER, ROUTES.HOMEPAGE);
  if (!isFilterEnabled) {
    return <></>;
  }
  // Hide datasets behind feature flag - end

  return (
    <>
      <Head>
        <title>cellxgene | Collections</title>
        {/* eslint-disable-next-line @next/next/no-page-custom-font -- required for font specs per mocks, revisit with #1685. */}
        <link
          href="https://fonts.googleapis.com/css?family=Roboto:400,500,700&amp;display=swap"
          rel="stylesheet"
        />
      </Head>
      <SideBar label="Filters">
        <Categories {...filterInstance} />
      </SideBar>
      <View>
        <FilteredCollectionsGrid tableInstance={tableInstance} />
      </View>
    </>
  );
}

/**
 * Create function that flattens and de-dupes array category values (in a single category/column) from dataset rows
 * grouped by collection. Used when aggregating and displaying dataset category values that are aggregated at the
 * collection level.
 * @returns Function that aggregates values across rows.
 */
function aggregateFn(): AggregatorFn<FilterableDataset> {
  // @param rows - array containing all values in a single category/column for all datasets grouped by collection
  return (_columnValues: CellValue[], rows: Array<Row<FilterableDataset>>) => {
    return [...new Set(rows.flat())];
  };
}

/**
 * Table cell component displaying multi-value cell values.
 * @param props - Cell-specific properties supplied from react-table.
 * @returns Array of DOM elements, one for each value in multi-value cell.
 */
function Cell(props: CellProps<FilterableDataset, string[]>): JSX.Element[] {
  // TODO(cc) useMemo?
  // TODO(cc) reuse with datasets
  const {
    cell: { value },
  } = props;
  return value.map((v: string) => <div key={v}>{v}</div>);
}

/**
 * Create function to be used by column.accessor in react-table column definition, for columns containing non-ontology
 * metadata (string) values. Used by react-table to determine filter and display values of column for each row
 * (in this case, the category value itself is filtered and displayed).
 * @param key - Object key of value to display in cell.
 * @returns Function that returns the values with the given key.
 */
function aggregatedCellAccessorFn(key: NonOntologyCategoryKey) {
  return (dataset: FilterableDataset) =>
    // @ts-expect-error -- TODO(cc) update filterableCollection to be required
    dataset.filterableCollection[`${key}Aggregated`];
}

/**
 * Create function to be used by column.accessor in react-table column definition, for columns containing ontology
 * metadata (ontology label and key) values. Used by react-table to determine filter and display values of column
 * for each row (in this case, the ontology label is filtered and displayed, not the full ontology object).
 // * @param key - Object key of value to display in cell.
 // * @returns Function that returns the array of ontology labels with the given key.
 */
function aggregatedOntologyCellAccessorFn(key: OntologyCategoryKey) {
  return (dataset: FilterableDataset): string[] => {
    // @ts-expect-error -- TODO(cc) update filterableCollection to be required
    return dataset.filterableCollection?.[`${key}Aggregated`].map(
      (o: Ontology) => o.label
    );
  };
}

/**
 * Group filterable datasets by collection.
 * @param filterableDatasets - Array of filterable datasets to group by their collection ID.
 * @returns Map of filterable datasets key by collection ID.
 */
function groupDatasetsByCollection(
  filterableDatasets: FilterableDataset[]
): Map<string, FilterableDataset[]> {
  return filterableDatasets.reduce(
    (accum: Map<string, FilterableDataset[]>, filterableDataset) => {
      const datasetsByCollectionId = accum.get(filterableDataset.collection_id);
      if (datasetsByCollectionId) {
        datasetsByCollectionId.push(filterableDataset);
      } else {
        accum.set(filterableDataset.collection_id, [filterableDataset]);
      }
      return accum;
    },
    new Map<string, FilterableDataset[]>()
  );
}

/**
 * Create filterable collections from aggregated dataset category values and add to each dataset in collection.
 * @param filterableDatasets - Filterable datasets to create filterable collections from.
 * @returns Array of updated filterable datasets, containing pointers to their corresponding filterable collections.
 */
function prepareData(
  filterableDatasets: FilterableDataset[]
): FilterableDataset[] {
  // Group datasets by collection to facilitate aggregation of dataset metadata per collection.
  const datasetsByCollectionId = groupDatasetsByCollection(filterableDatasets);

  // Aggregate metadata for each collection and update on each dataset.
  const groupedFilterableDatasets = [...datasetsByCollectionId.values()].map(
    (filterableDatasets: FilterableDataset[]) => {
      // Create model of collection category values by aggregating the values in each category of each dataset in
      // collection.
      const filterableCollection =
        createFilterableCollection(filterableDatasets);

      // Add aggregated collection category values to each dataset
      return filterableDatasets.map((filterableDataset: FilterableDataset) => ({
        ...filterableDataset,
        filterableCollection,
      }));
    }
  );

  // Flatten the array of filterable datasets array.
  return groupedFilterableDatasets.flat();
}

/**
 * Create model of collection category values by aggregating the values in each category of each dataset in collection.
 * @param filterableDatasets - Datasets in the collection to aggregate category values over.
 * @returns Collection object containing aggregated category values from given filterable datasets.
 */
function createFilterableCollection(
  filterableDatasets: FilterableDataset[]
): FilterableCollection {
  return Object.values(CATEGORY_KEY).reduce(
    (accum: FilterableCollection, categoryKey: CategoryKey) => {
      // @ts-expect-error -- TODO(cc) revisit
      accum[`${categoryKey}Aggregated`] = aggregateCategory(
        categoryKey,
        filterableDatasets
      );
      return accum;
    },
    {} as FilterableCollection
  );
}

/**
 * Determine the set of category values that exist in the given datasets for the given category.
 * @param categoryKey - Key of the category to aggregate values of.
 * @param filterableDatasets - Datasets to aggregate category values of.
 * @returns Array of aggregated category values for the given category.
 */
function aggregateCategory(
  categoryKey: CategoryKey,
  filterableDatasets: FilterableDataset[]
): (Ontology | IS_PRIMARY_DATA)[] {
  // TODO(cc) return type
  const categoryValuesSet = filterableDatasets.reduce(
    (
      accum: Set<Ontology | IS_PRIMARY_DATA>,
      filterableDataset: FilterableDataset
    ) => {
      const categoryValue = filterableDataset[categoryKey];
      if (!Array.isArray(categoryValue)) {
        accum.add(categoryValue);
      } else {
        categoryValue.forEach((categoryValue) => accum.add(categoryValue));
      }
      return accum;
    },
    new Set<Ontology | IS_PRIMARY_DATA>()
  );
  return [...categoryValuesSet];
}
