import Head from "next/head";
import React, { useMemo } from "react";
import { Column, useFilters, useSortBy, useTable } from "react-table";
import { PluralizedMetadataLabel } from "src/common/constants/metadata";
import { ROUTES } from "src/common/constants/routes";
import { FEATURES } from "src/common/featureFlags/features";
import { useCategoryFilter } from "src/common/hooks/useCategoryFilter";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { useFetchCollectionRows } from "src/common/queries/filter";
import { CollectionsGrid } from "src/components/Collections/components/Grid/components/CollectionsGrid/style";
import Filter from "src/components/common/Filter";
import {
  CATEGORY_KEY,
  CellPropsValue,
  CollectionRow,
  RowPropsValue,
} from "src/components/common/Filter/common/entities";
import { ontologyCellAccessorFn } from "src/components/common/Filter/common/utils";
import { GridHero } from "src/components/common/Grid/components/Hero";
import LinkCell from "src/components/common/Grid/components/LinkCell";
import NTagCell from "src/components/common/Grid/components/NTagCell";
import SideBar from "src/components/common/SideBar";
import { View } from "src/views/globalStyle";

// Collection ID object key
const COLLECTION_ID = "id";

// Collection name object key
const COLLECTION_NAME = "name";

// Key identifying recency sort by column
const COLUMN_ID_RECENCY = "recency";

export default function Collections(): JSX.Element {
  // Filterable collection datasets joined from datasets index and collections index responses.
  const { error, loading, rows: collectionRows } = useFetchCollectionRows();

  // Column configuration backing table.
  const columnConfig: Column<CollectionRow>[] = useMemo(
    () => [
      {
        Cell: ({ row }: RowPropsValue) => {
          return (
            <LinkCell
              data-test-id="collection-link"
              url={ROUTES.COLLECTION.replace(":id", row.values.id)}
              value={row.values.name}
            />
          );
        },
        Header: "Collection",
        accessor: COLLECTION_NAME,
      },
      {
        Cell: ({ value }: CellPropsValue) => (
          <NTagCell label={PluralizedMetadataLabel.TISSUE} values={value} />
        ),
        Header: "Tissue",
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.TISSUE),
        filter: "includesSome",
        id: CATEGORY_KEY.TISSUE,
      },
      {
        Cell: ({ value }: CellPropsValue) => (
          <NTagCell label={PluralizedMetadataLabel.DISEASE} values={value} />
        ),
        Header: "Disease",
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.DISEASE),
        filter: "includesSome",
        id: CATEGORY_KEY.DISEASE,
      },
      {
        Cell: ({ value }: CellPropsValue) => (
          <NTagCell label={PluralizedMetadataLabel.ORGANISM} values={value} />
        ),
        Header: "Organism",
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.ORGANISM),
        filter: "includesSome",
        id: CATEGORY_KEY.ORGANISM,
      },
      // Hidden, required for sorting
      {
        // Sort by revised_at if specified otherwise published_at.
        accessor: (dataset: CollectionRow): number =>
          dataset.revised_at ?? dataset.published_at,
        id: COLUMN_ID_RECENCY,
      },
      // Hidden, required for accessing collection ID via row.values, for building link to collection detail page.
      {
        accessor: COLLECTION_ID,
      },
      // Hidden, required for filter.
      {
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.ASSAY),
        filter: "includesSome",
        id: CATEGORY_KEY.ASSAY,
      },
      // Hidden, required for filter.
      {
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.CELL_TYPE),
        filter: "includesSome",
        id: CATEGORY_KEY.CELL_TYPE,
      },
      // Hidden, required for filter.
      {
        accessor: (dataset: CollectionRow) => dataset.is_primary_data,
        filter: "includesSome",
        id: CATEGORY_KEY.IS_PRIMARY_DATA,
      },
      // Hidden, required for filter.
      {
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.SEX),
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
      data: collectionRows,
      initialState: {
        // Only display tissue, disease and organism values.
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
    useSortBy
  );

  // Filter init.
  const {
    preFilteredRows,
    rows,
    setFilter,
    state: { filters },
  } = tableInstance;
  const filterInstance = useCategoryFilter(preFilteredRows, filters, setFilter);

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
      {error || loading ? null : (
        <>
          <SideBar label="Filters">
            <Filter {...filterInstance} />
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
      )}
    </>
  );
}
