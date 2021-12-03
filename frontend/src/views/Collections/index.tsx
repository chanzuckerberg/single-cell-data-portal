import Head from "next/head";
import React, { useMemo, useState } from "react";
import { Column, useFilters, useSortBy, useTable } from "react-table";
import { PluralizedMetadataLabel } from "src/common/constants/metadata";
import { ROUTES } from "src/common/constants/routes";
import { FEATURES } from "src/common/featureFlags/features";
import { useCategoryFilter } from "src/common/hooks/useCategoryFilter";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { fetchCollectionRows } from "src/common/queries/filter";
import Categories from "src/components/Categories";
import { CollectionsGrid } from "src/components/Collections/components/Grid/components/CollectionsGrid/style";
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
const COLLECTION_ID = "collection_id";

// Collection name object key
const COLLECTION_NAME = "name";

// Key identifying recency sort by column
const COLUMN_ID_RECENCY = "recency";

export default function Collections(): JSX.Element {
  // Filterable collection datasets joined from datasets index and collections index responses.
  const [filterableCollectionDatasets] = useState<CollectionRow[]>(
    fetchCollectionRows()
  );

  // Column configuration backing table.
  const columnConfig: Column<CollectionRow>[] = useMemo(
    // TODO(cc) remove cell and header from hidden cols (same for datasets)
    // TODO(cc) share ontology accessor with datasets
    () => [
      // Hidden, required for sorting
      {
        // Sort by revised_at if specified otherwise published_at.
        accessor: (dataset: CollectionRow): number =>
          dataset.revised_at ?? dataset.published_at,
        id: COLUMN_ID_RECENCY,
      },
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
      {
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.ASSAY),
        filter: "includesSome",
        id: CATEGORY_KEY.ASSAY,
      },
      {
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.CELL_TYPE),
        filter: "includesSome",
        id: CATEGORY_KEY.CELL_TYPE,
      },
      {
        accessor: (dataset: CollectionRow) => dataset.is_primary_data,
        filter: "includesSome",
        id: CATEGORY_KEY.IS_PRIMARY_DATA,
      },
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
      data: filterableCollectionDatasets,
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
