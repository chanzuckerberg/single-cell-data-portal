import Head from "next/head";
import React, { FC, useMemo } from "react";
import {
  AggregatorFn,
  CellProps,
  CellValue,
  Column,
  Row,
  useFilters,
  useGroupBy,
  useTable,
} from "react-table";
import { ROUTES } from "src/common/constants/routes";
import {
  FilterableCollection,
  FilterableDataset,
  IS_PRIMARY_DATA,
  Ontology,
} from "src/common/entities";
import { FEATURES } from "src/common/featureFlags/features";
import {
  CategoryKey,
  CATEGORY_KEY,
  OntologyCategoryKey,
  useFacetedFilter,
} from "src/common/hooks/useFacetedFilter";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import Categories from "src/components/Categories";
import FilteredCollectionsGrid from "src/components/Collections/components/Grid/components/FilteredCollectionsGrid";
import SideBar from "src/components/common/SideBar";
import { View } from "src/views/globalStyle";
import filterableDatasets from "../../../tests/features/fixtures/datasets/filterable-datasets";

const Collections: FC = () => {
  // Column configuration backing table.
  // TODO(cc) function for creating aggregate headers and function for creating filtering headers
  const columnConfig: Column<FilterableDataset>[] = useMemo(
    () => [
      // Required for grouping datasets by collections, not displayed.
      {
        Header: "Collection ID",
        accessor: "collection_id",
      },
      // Collection name, aggregated across datasets.
      {
        Header: "Collection",
        accessor: "collection_name", // TODO(cc) use constant here?
        aggregate: aggregateFn(),
      },
      // // Tissue for specific to each dataset row. Used for filtering datasets by the selected category values. Not
      // // displayed.
      // {
      //   Header: "Tissue",
      //   accessor: ontologyCellAccessorFn(CATEGORY_KEY.TISSUE),
      //   filter: "includesSome",
      //   id: CATEGORY_KEY.TISSUE,
      // },
      // Unfiltered, aggregated tissue values across datasets in a collection. Used for displaying tissue values for
      // an aggregated "materialized" collection row.
      {
        Cell: Cell,
        Header: "Tissue",
        accessor: filteredCollectionOntologyAccessorFn(CATEGORY_KEY.TISSUE),
        aggregate: aggregateFn(), // TODO(cc) can this be made into some kind of typed or named fn? same with accessor.
        filter: "includesSome",
        id: `${CATEGORY_KEY.TISSUE}`,
      },
      // Disease for specific to each dataset row. Used for filtering datasets by the selected category values. Not
      // displayed.
      // {
      //   Header: "Disease",
      //   accessor: ontologyCellAccessorFn(CATEGORY_KEY.DISEASE),
      //   filter: "includesSome",
      //   id: `${CATEGORY_KEY.DISEASE}`,
      // },
      // Unfiltered, aggregated disease values across datasets in a collection. Used for displaying disease values for
      // an aggregated "materialized" collection row.
      {
        Cell: Cell,
        Header: "Disease",
        accessor: filteredCollectionOntologyAccessorFn(CATEGORY_KEY.DISEASE),
        aggregate: aggregateFn(),
        filter: "includesSome",
        id: `${CATEGORY_KEY.DISEASE}`, // TODO(cc) function for creating aggregate key name
      },
      // Assay for specific to each dataset row. Used for filtering datasets by the selected category values. Not
      // displayed.
      // {
      //   Header: "Assay",
      //   accessor: ontologyCellAccessorFn(CATEGORY_KEY.ASSAY),
      //   filter: "includesSome",
      //   id: CATEGORY_KEY.ASSAY,
      // },
      // Unfiltered, aggregated disease values across datasets in a collection. Used for displaying disease values for
      // an aggregated "materialized" collection row.
      {
        Cell: Cell,
        Header: "Assay",
        accessor: filteredCollectionOntologyAccessorFn(CATEGORY_KEY.ASSAY),
        aggregate: aggregateFn(),
        filter: "includesSome",
        id: `${CATEGORY_KEY.ASSAY}`,
      },
      // Organism for specific to each dataset row. Used for filtering datasets by the selected category values. Not
      // displayed.
      // {
      //   Header: "Organism",
      //   accessor: ontologyCellAccessorFn(CATEGORY_KEY.ORGANISM),
      //   filter: "includesSome",
      //   id: CATEGORY_KEY.ORGANISM,
      // },
      // Unfiltered, aggregated organism values across datasets in a collection. Used for displaying organism values for
      // an aggregated "materialized" collection row.
      {
        Cell: Cell,
        Header: "Organism",
        accessor: filteredCollectionOntologyAccessorFn(CATEGORY_KEY.ORGANISM),
        aggregate: aggregateFn(),
        filter: "includesSome",
        id: `${CATEGORY_KEY.ORGANISM}`,
      },
      // Header spec for filtering and aggregating cell type.
      // {
      //   Header: "Cell Type",
      //   accessor: ontologyCellAccessorFn(CATEGORY_KEY.CELL_TYPE),
      //   filter: "includesSome",
      //   id: CATEGORY_KEY.CELL_TYPE,
      // },
      // Unfiltered, aggregated cell type values across datasets in a collection. Used for displaying organism values
      // for an aggregated "materialized" collection row.
      {
        Cell: Cell,
        Header: "Cell Type",
        accessor: filteredCollectionOntologyAccessorFn(CATEGORY_KEY.CELL_TYPE),
        aggregate: aggregateFn(),
        filter: "includesSome",
        id: `${CATEGORY_KEY.CELL_TYPE}`,
      },
      // Header spec for filtering and aggregating primary data.
      // {
      //   Header: "Primary Data",
      //   accessor: CATEGORY_KEY.IS_PRIMARY_DATA,
      //   filter: "includesSome",
      //   id: CATEGORY_KEY.IS_PRIMARY_DATA,
      // },
      // Unfiltered, aggregated primary data values across datasets in a collection. Used for displaying organism values
      // for an aggregated "materialized" collection row.
      {
        Cell: Cell,
        Header: "Primary Data",
        accessor: (dataset: FilterableDataset) =>
          dataset.filterableCollection?.[
            `${CATEGORY_KEY.IS_PRIMARY_DATA}Aggregated`
          ], // TODO(cc) update ontology accessor to use is-a on ontology or string and return value from there
        aggregate: aggregateFn(),
        filter: "includesSome",
        id: `${CATEGORY_KEY.IS_PRIMARY_DATA}`,
      },
      // Header spec for filtering and aggregating sex.
      // {
      //   Header: "Sex",
      //   accessor: ontologyCellAccessorFn(CATEGORY_KEY.SEX),
      //   filter: "includesSome",
      //   id: CATEGORY_KEY.SEX,
      // },
      // Unfiltered, aggregated sex values across datasets in a collection. Used for displaying organism values
      // for an aggregated "materialized" collection row.
      {
        Cell: Cell,
        Header: "Sex",
        accessor: filteredCollectionOntologyAccessorFn(CATEGORY_KEY.SEX),
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
        groupBy: ["collection_id"],
        // Only display aggregated tissue, disease and organism values.
        hiddenColumns: [
          "collection_id", // TODO(cc) constant?
          // CATEGORY_KEY.ASSAY,
          // `${CATEGORY_KEY.ORGANISM}Aggregated`, // TODO(cc) remove for go-live (keep for manual testing purposes)
          // CATEGORY_KEY.DISEASE,
          // CATEGORY_KEY.CELL_TYPE,
          // `${CATEGORY_KEY.CELL_TYPE}Aggregated`, // TODO(cc) remove for go-live (keep for manual testing purposes)
          // CATEGORY_KEY.IS_PRIMARY_DATA,
          // `${CATEGORY_KEY.IS_PRIMARY_DATA}Aggregated`, // TODO(cc) remove for go-live (keep for manual testing purposes)
          // CATEGORY_KEY.ORGANISM,
          // `${CATEGORY_KEY.ORGANISM}Aggregated`, // TODO(cc) remove for go-live (keep for manual testing purposes)
          // CATEGORY_KEY.SEX,
          // `${CATEGORY_KEY.SEX}Aggregated`, // TODO(cc) remove for go-live (keep for manual testing purposes)
          // CATEGORY_KEY.TISSUE,
        ],
      },
    },
    useFilters,
    useGroupBy
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
    "collection_id" // TODO(cc) constant
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
};

/**
 * Create function that de-dupes aggregated values into single array. Used when displaying dataset metadata values that
 * are aggregated at the collection level.
 * @returns Function that aggregates values across rows.
 */
function aggregateFn(): AggregatorFn<FilterableDataset> {
  // TODO(cc) review types
  return (_columnValues: CellValue[], rows: Array<Row<FilterableDataset>>) => [
    ...new Set(rows.flat()),
  ];
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
 * TODO(cc) return type etc
 */
function filteredCollectionOntologyAccessorFn(key: OntologyCategoryKey) {
  return (dataset: FilterableDataset) => {
    // TODO(cc) ? optional filterableCollection - maybe filterable collection should be its own type and collection should be on the side?
    return dataset.filterableCollection?.[`${key}Aggregated`].map(
      (o: Ontology) => o.label
    );
  };
}

/**
 * TODO(cc)
 * @param filterableDatasets
 */
function groupDatasetsByCollection( // TODO(cc) - duped
  filterableDatasets: FilterableDataset[]
): Map<string, FilterableDataset[]> {
  return filterableDatasets.reduce(
    (accum: Map<string, FilterableDataset[]>, filterableDataset) => {
      let datasetsByCollectionId = accum.get(filterableDataset.collection_id);
      if (!datasetsByCollectionId) {
        datasetsByCollectionId = [];
        accum.set(filterableDataset.collection_id, datasetsByCollectionId);
      }
      datasetsByCollectionId.push(filterableDataset);
      return accum;
    },
    new Map<string, FilterableDataset[]>()
  );
}

/**
 * Create function to be used by column.accessor in react-table column definition, for columns containing ontology
 * metadata (ontology label and key) values.
 * @param key - Object key of value to display in cell.
 * @returns Function that returns the value with the given key, to display in a cell.
 */
// function ontologyCellAccessorFn( TODO(cc) remove
//   key: OntologyCategoryKey
// ): (dataset: FilterableDataset) => string[] {
//   // TODO(cc) reuse with datasets
//   return (dataset: FilterableDataset) =>
//     dataset[key].map((o: Ontology) => o.label);
// }

/**
 * TODO(cc)
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
 * @param filterableDatasets - Datasets in the collection to aggregate metadata of.
 * @returns Collection object containing aggregated category values across. TODO(cc)
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

export default Collections;
