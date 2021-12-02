import Head from "next/head";
import React, { useMemo, useState } from "react";
import {
  AggregatorFn,
  CellValue,
  Column,
  Row,
  useFilters,
  useGroupBy,
  useSortBy,
  useTable,
} from "react-table";
import { PluralizedMetadataLabel } from "src/common/constants/metadata";
import { ROUTES } from "src/common/constants/routes";
import { Ontology } from "src/common/entities";
import { FEATURES } from "src/common/featureFlags/features";
import {
  NonOntologyCategoryKey,
  OntologyCategoryKey,
  useCategoryFilter,
} from "src/common/hooks/useCategoryFilter";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { fetchCollectionRows } from "src/common/queries/filterable-datasets";
import Categories from "src/components/Categories";
import { CollectionsGrid } from "src/components/Collections/components/Grid/components/CollectionsGrid/style";
import {
  CATEGORY_KEY,
  CellPropsValue,
  CollectionRow,
  RowPropsValue,
} from "src/components/common/Filter/common/entities";
import Cell from "src/components/common/Grid/components/Cell";
import { GridHero } from "src/components/common/Grid/components/Hero";
import LinkCell from "src/components/common/Grid/components/LinkCell";
import NTagCell from "src/components/common/Grid/components/NTagCell";
import SideBar from "src/components/common/SideBar";
import { View } from "src/views/globalStyle";

// Collection ID object key
const COLLECTION_ID = "collection_id";

// Collection name object key
const COLLECTION_NAME = "collection_name";

// Key identifying recency sort by column
const COLUMN_ID_RECENCY = "recency";

export default function Collections(): JSX.Element {
  // Filterable collection datasets joined from datasets index and collections index responses.
  const [filterableCollectionDatasets] = useState<CollectionRow[]>(
    fetchCollectionRows()
  );

  // Column configuration backing table.
  const columnConfig: Column<CollectionRow>[] = useMemo(
    () => [
      // Hidden, required for grouping datasets by collections
      {
        accessor: COLLECTION_ID,
      },
      // Hidden, required for sorting TODO(cc) this should be for materialized collection row
      {
        // Sort by revised_at if specified otherwise published_at.
        accessor: (dataset: CollectionRow): number =>
          dataset.revised_at ?? dataset.published_at,
        id: COLUMN_ID_RECENCY,
      },
      // Collection name, aggregated across datasets to roll up into single value for each collection ID group.
      {
        Cell: ({ row }: RowPropsValue) => {
          return (
            <LinkCell
              data-test-id="collection-link"
              url={ROUTES.COLLECTION.replace(":id", row.values.collection_id)}
              value={row.values.collection_name[0]}
            />
          );
        },
        Header: "Collection",
        accessor: COLLECTION_NAME,
        aggregate: aggregateFn(),
      },
      {
        Cell: ({ value }: CellPropsValue) => (
          <NTagCell label={PluralizedMetadataLabel.TISSUE} values={value} />
        ),
        Header: "Tissue",
        accessor: aggregatedOntologyCellAccessorFn(CATEGORY_KEY.TISSUE),
        aggregate: aggregateFn(),
        filter: "includesSome",
        id: CATEGORY_KEY.TISSUE,
      },
      {
        Cell: ({ value }: CellPropsValue) => (
          <NTagCell label={PluralizedMetadataLabel.DISEASE} values={value} />
        ),
        Header: "Disease",
        accessor: aggregatedOntologyCellAccessorFn(CATEGORY_KEY.DISEASE),
        aggregate: aggregateFn(),
        filter: "includesSome",
        id: CATEGORY_KEY.DISEASE,
      },
      {
        Cell: ({ value }: CellPropsValue) => (
          <NTagCell label={PluralizedMetadataLabel.ASSAY} values={value} />
        ),
        Header: "Assay",
        accessor: aggregatedOntologyCellAccessorFn(CATEGORY_KEY.ASSAY),
        aggregate: aggregateFn(),
        filter: "includesSome",
        id: CATEGORY_KEY.ASSAY,
      },
      {
        Cell: ({ value }: CellPropsValue) => (
          <NTagCell label={PluralizedMetadataLabel.ORGANISM} values={value} />
        ),
        Header: "Organism",
        accessor: aggregatedOntologyCellAccessorFn(CATEGORY_KEY.ORGANISM),
        aggregate: aggregateFn(),
        filter: "includesSome",
        id: CATEGORY_KEY.ORGANISM, // TODO(cc) check for `` in datasets
      },
      {
        Cell: ({ value }: CellPropsValue) => (
          <NTagCell label={PluralizedMetadataLabel.CELL_TYPE} values={value} />
        ), // TODO(cc) remove cell and header from hidden cols (same for datasets)
        Header: "Cell Type",
        accessor: aggregatedOntologyCellAccessorFn(CATEGORY_KEY.CELL_TYPE),
        aggregate: aggregateFn(),
        filter: "includesSome",
        id: CATEGORY_KEY.CELL_TYPE,
      },
      {
        Cell: Cell,
        Header: "Primary Data",
        accessor: aggregatedCellAccessorFn(CATEGORY_KEY.IS_PRIMARY_DATA),
        aggregate: aggregateFn(),
        filter: "includesSome",
        id: CATEGORY_KEY.IS_PRIMARY_DATA,
      },
      {
        Cell: Cell,
        Header: "Sex",
        accessor: aggregatedOntologyCellAccessorFn(CATEGORY_KEY.SEX),
        aggregate: aggregateFn(),
        filter: "includesSome",
        id: CATEGORY_KEY.SEX,
      },
    ],
    []
  );

  // Table init
  const tableInstance = useTable<CollectionRow>(
    {
      columns: columnConfig,
      data: filterableCollectionDatasets,
      initialState: {
        groupBy: [COLLECTION_ID],
        // Only display aggregated tissue, disease and organism values.
        hiddenColumns: [
          COLLECTION_ID,
          COLUMN_ID_RECENCY,
          CATEGORY_KEY.ASSAY,
          CATEGORY_KEY.CELL_TYPE,
          CATEGORY_KEY.IS_PRIMARY_DATA,
          CATEGORY_KEY.SEX,
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
    rows,
    setFilter,
    state: { filters },
  } = tableInstance;
  const filterInstance = useCategoryFilter(
    // @ts-expect-error -- TODO(cc) revisit
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
        {!rows || rows.length === 0 ? (
          <GridHero>
            <h3>No Results</h3>
            <p>There are no collections matching those filters.</p>
          </GridHero>
        ) : (
          <CollectionsGrid tableInstance={tableInstance} />
        )}
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
function aggregateFn(): AggregatorFn<CollectionRow> {
  // TODO(cc) can this just be the outside function?
  // @param rows - array containing all values in a single category/column for all datasets grouped by collection
  return (_columnValues: CellValue[], rows: Array<Row<CollectionRow>>) => {
    return [...new Set(rows.flat())];
  };
}

/**
 * Create function to be used by column.accessor in react-table column definition, for columns containing non-ontology
 * metadata (string) values. Used by react-table to determine filter and display values of column for each row
 * (in this case, the category value itself is filtered and displayed).
 * @param key - Object key of value to display in cell.
 * @returns Function that returns the values with the given key.
 */
function aggregatedCellAccessorFn(key: NonOntologyCategoryKey) {
  return (dataset: CollectionRow) => dataset[`${key}Aggregated`];
}

/**
 * Create function to be used by column.accessor in react-table column definition, for columns containing ontology
 * metadata (ontology label and key) values. Used by react-table to determine filter and display values of column
 * for each row (in this case, the ontology label is filtered and displayed, not the full ontology object).
 // * @param key - Object key of value to display in cell.
 // * @returns Function that returns the array of ontology labels with the given key.
 */
function aggregatedOntologyCellAccessorFn(key: OntologyCategoryKey) {
  return (dataset: CollectionRow): string[] =>
    dataset[`${key}Aggregated`].map((o: Ontology) => o.label);
}
