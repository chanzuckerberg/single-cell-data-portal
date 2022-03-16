import loadable from "@loadable/component";
import Head from "next/head";
import React, { useMemo } from "react";
import { Column, useFilters, useSortBy, useTable } from "react-table";
import { PLURALIZED_METADATA_LABEL } from "src/common/constants/metadata";
import { ROUTES } from "src/common/constants/routes";
import { FEATURES } from "src/common/featureFlags/features";
import {
  CategoryKey,
  useCategoryFilter,
} from "src/common/hooks/useCategoryFilter";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { useFetchCollectionRows } from "src/common/queries/filter";
import { useExplainTombstoned } from "src/components/Collections/common/utils";
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
import { Title } from "src/components/common/Grid/components/Title";
import SideBar from "src/components/common/SideBar";
import { View } from "src/views/globalStyle";

/**
 * Collection ID object key.
 */
const COLLECTION_ID = "id";

/**
 * Collection name object key.
 */
const COLLECTION_NAME = "name";

/**
 * Collection summary citation object key.
 */
const COLLECTION_SUMMARY_CITATION = "summaryCitation";

/**
 * Key identifying recency sort by column.
 */
const COLUMN_ID_RECENCY = "recency";

/**
 * Recency object key.
 */
const RECENCY = "recency";

/**
 * Gene sets CSV upload functionality, available if gene sets feature flag is enabled
 */
const AsyncUploadCSV = loadable(
  () =>
    /*webpackChunkName: 'src/components/UploadCSV' */ import(
      "src/components/UploadCSV"
    )
);

export default function Collections(): JSX.Element {
  // Pop toast if user has been redirected from a tombstoned collection.
  useExplainTombstoned();

  // Determine if gene sets functionality is available to user.
  const isGeneSetsOn = useFeatureFlag(FEATURES.GENE_SETS);

  // Filterable collection datasets joined from datasets index and collections index responses.
  const { isError, isLoading, rows: collectionRows } = useFetchCollectionRows();

  // Column configuration backing table.
  const columnConfig: Column<CollectionRow>[] = useMemo(
    () => [
      {
        Cell: ({ row }: RowPropsValue<CollectionRow>) => {
          return (
            <Title>
              <LinkCell
                data-test-id="collection-link"
                url={ROUTES.COLLECTION.replace(":id", row.values.id)}
                value={row.values.name}
              />
            </Title>
          );
        },
        Header: "Collection",
        accessor: COLLECTION_NAME,
      },
      {
        Cell: ({ row }: RowPropsValue<CollectionRow>) => {
          return <div>{row.values.summaryCitation || "No publication"}</div>;
        },
        Header: "Publication",
        accessor: COLLECTION_SUMMARY_CITATION,
      },
      {
        Cell: ({ value }: CellPropsValue<string[]>) => (
          <NTagCell label={PLURALIZED_METADATA_LABEL.TISSUE} values={value} />
        ),
        Header: "Tissue",
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.TISSUE),
        filter: "includesSome",
        id: CATEGORY_KEY.TISSUE,
      },
      {
        Cell: ({ value }: CellPropsValue<string[]>) => (
          <NTagCell label={PLURALIZED_METADATA_LABEL.DISEASE} values={value} />
        ),
        Header: "Disease",
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.DISEASE),
        filter: "includesSome",
        id: CATEGORY_KEY.DISEASE,
      },
      {
        Cell: ({ value }: CellPropsValue<string[]>) => (
          <NTagCell label={PLURALIZED_METADATA_LABEL.ORGANISM} values={value} />
        ),
        Header: "Organism",
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.ORGANISM),
        filter: "includesSome",
        id: CATEGORY_KEY.ORGANISM,
      },
      // Hidden, required for sorting
      {
        accessor: RECENCY,
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
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.ETHNICITY),
        filter: "includesSome",
        id: CATEGORY_KEY.ETHNICITY,
      },
      // Hidden, required for filter.
      {
        accessor: (collectionRow: CollectionRow) =>
          collectionRow.is_primary_data,
        filter: "includesSome",
        id: CATEGORY_KEY.IS_PRIMARY_DATA,
      },
      // Hidden, required for filter.
      {
        accessor: CATEGORY_KEY.PUBLICATION_AUTHORS,
        filter: "includesSome",
        id: CATEGORY_KEY.PUBLICATION_AUTHORS,
      },
      // Hidden, required for filter.
      {
        accessor: CATEGORY_KEY.PUBLICATION_DATE_VALUES,
        filter: "includesSome",
        id: CATEGORY_KEY.PUBLICATION_DATE_VALUES,
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
          CATEGORY_KEY.ETHNICITY,
          CATEGORY_KEY.IS_PRIMARY_DATA,
          CATEGORY_KEY.PUBLICATION_AUTHORS,
          CATEGORY_KEY.PUBLICATION_DATE_VALUES,
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

  // Determine the set of categories to display for the datasets view.
  const categories = useMemo<Set<CATEGORY_KEY>>(() => {
    return Object.values(CATEGORY_KEY)
      .filter(
        (categoryKey: CategoryKey) =>
          categoryKey !== CATEGORY_KEY.CELL_COUNT &&
          categoryKey !== CATEGORY_KEY.MEAN_GENES_PER_CELL
      )
      .reduce((accum, categoryKey: CategoryKey) => {
        accum.add(categoryKey);
        return accum;
      }, new Set<CATEGORY_KEY>());
  }, []);

  // Filter init.
  const {
    preFilteredRows,
    rows,
    setFilter,
    state: { filters },
  } = tableInstance;
  const filterInstance = useCategoryFilter(
    preFilteredRows,
    categories,
    filters,
    setFilter
  );

  // Hide datasets behind feature flag - start
  const isFilterEnabled = useFeatureFlag(FEATURES.FILTER);
  if (!isFilterEnabled) {
    return <></>;
  }
  // Hide datasets behind feature flag - end

  return (
    <>
      <Head>
        <title>cellxgene | Collections</title>
      </Head>
      {isError || isLoading ? null : (
        <>
          <SideBar label="Filters" isOpen>
            <Filter {...filterInstance} />
          </SideBar>
          <View>
            {
              // (thuang): TEMP. Remove when we do https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/corpora-data-portal/917
              isGeneSetsOn && <AsyncUploadCSV />
            }
            {!rows || rows.length === 0 ? (
              <GridHero>
                <h3>No Results</h3>
                <p>There are no collections matching those filters.</p>
              </GridHero>
            ) : (
              // @ts-expect-error -- revisit tableInstance typing
              <CollectionsGrid tableInstance={tableInstance} />
            )}
          </View>
        </>
      )}
    </>
  );
}
