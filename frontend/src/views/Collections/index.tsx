import Head from "next/head";
import React, { useEffect, useMemo, useState } from "react";
import { Column, Filters, useFilters, useSortBy, useTable } from "react-table";
import { PLURALIZED_METADATA_LABEL } from "src/common/constants/metadata";
import { ROUTES } from "src/common/constants/routes";
import { FEATURES } from "src/common/featureFlags/features";
import { useCategoryFilter } from "src/common/hooks/useCategoryFilter";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { useSessionStorage } from "src/common/hooks/useSessionStorage";
import { useFetchCollectionRows } from "src/common/queries/filter";
import { KEYS } from "src/common/sessionStorage/set";
import { useExplainTombstoned } from "src/components/Collections/common/utils";
import { CollectionsGrid } from "src/components/Collections/components/Grid/components/CollectionsGrid/style";
import Filter from "src/components/common/Filter";
import {
  CATEGORY_FILTER_ID,
  CellPropsValue,
  CollectionRow,
  RowPropsValue,
} from "src/components/common/Filter/common/entities";
import { ontologyLabelCellAccessorFn } from "src/components/common/Filter/common/utils";
import DiseaseCell from "src/components/common/Grid/components/DiseaseCell";
import HeaderCell from "src/components/common/Grid/components/HeaderCell";
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

export default function Collections(): JSX.Element {
  // Pop toast if user has been redirected from a tombstoned collection.
  useExplainTombstoned();

  // Filterable collection datasets joined from datasets index and collections index responses.
  const { isError, isLoading, rows: collectionRows } = useFetchCollectionRows();
  const [rowCount, setRowCount] = useState<number>();
  const totalCount = collectionRows.length;

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
        Header: (
          <HeaderCell
            label={"Collections"}
            rowCount={rowCount}
            totalCount={totalCount}
          />
        ),
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
        accessor: ontologyLabelCellAccessorFn("tissue"),
        filter: "includesSome", // TODO(cc) remove filter with 2569? but keep tissue column for display (tissue filter below will filter on term IDs but tissue column here will display labels)
        id: CATEGORY_FILTER_ID.TISSUE_DEPRECATED,
      },
      {
        Cell: ({ value }: CellPropsValue<string[]>) => (
          <DiseaseCell
            label={PLURALIZED_METADATA_LABEL.DISEASE}
            values={value}
          />
        ),
        Header: "Disease",
        accessor: ontologyLabelCellAccessorFn("disease"),
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.DISEASE,
      },
      {
        Cell: ({ value }: CellPropsValue<string[]>) => (
          <NTagCell label={PLURALIZED_METADATA_LABEL.ORGANISM} values={value} />
        ),
        Header: "Organism",
        accessor: ontologyLabelCellAccessorFn("organism"),
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.ORGANISM,
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
        accessor: ontologyLabelCellAccessorFn("assay"),
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.ASSAY,
      },
      // Hidden, required for filter.
      {
        accessor: ontologyLabelCellAccessorFn("cell_type"),
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.CELL_TYPE_DEPRECATED,
      },
      // Hidden, required for filter.
      {
        accessor: "cell_type_ancestors",
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.CELL_TYPE,
      },
      // Hidden, required for filter.
      {
        accessor: ontologyLabelCellAccessorFn("ethnicity"),
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.ETHNICITY,
      },
      {
        accessor: "development_stage_ancestors",
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
      },
      // Hidden, required for filter.
      {
        accessor: "publicationAuthors",
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.PUBLICATION_AUTHORS,
      },
      // Hidden, required for filter.
      {
        accessor: "publicationDateValues",
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.PUBLICATION_DATE_VALUES,
      },
      // Hidden, required for filter.
      {
        accessor: ontologyLabelCellAccessorFn("sex"),
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.SEX,
      },
      // Hidden, required for filter.
      {
        accessor: "tissueCalculated",
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.TISSUE_CALCULATED,
      },
    ],
    [rowCount, totalCount]
  );

  // Handle initial filter state and save of filter state beyond component scope.
  const [initialFilters, storeFilters] = useSessionStorage<
    Filters<CollectionRow>
  >(KEYS.FILTER_COLLECTIONS, []);

  // Table init
  const tableInstance = useTable<CollectionRow>(
    {
      columns: columnConfig,
      data: collectionRows,
      initialState: {
        filters: initialFilters,
        // Only display tissue, disease and organism values.
        hiddenColumns: [
          COLLECTION_ID,
          COLUMN_ID_RECENCY,
          CATEGORY_FILTER_ID.ASSAY,
          CATEGORY_FILTER_ID.CELL_TYPE_DEPRECATED,
          CATEGORY_FILTER_ID.CELL_TYPE,
          CATEGORY_FILTER_ID.ETHNICITY,
          CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
          CATEGORY_FILTER_ID.PUBLICATION_AUTHORS,
          CATEGORY_FILTER_ID.PUBLICATION_DATE_VALUES,
          CATEGORY_FILTER_ID.SEX,
          CATEGORY_FILTER_ID.TISSUE_CALCULATED,
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
  const isFilterEnabled = useFeatureFlag(FEATURES.FILTER); // TODO(cc) remove with #2569.
  const categories = useMemo<Set<CATEGORY_FILTER_ID>>(() => {
    return Object.values(CATEGORY_FILTER_ID)
      .filter((categoryFilterId: CATEGORY_FILTER_ID) => {
        if (
          categoryFilterId === CATEGORY_FILTER_ID.CELL_COUNT ||
          categoryFilterId == CATEGORY_FILTER_ID.GENE_COUNT
        ) {
          return false;
        }
        return !(
          (categoryFilterId === CATEGORY_FILTER_ID.TISSUE_CALCULATED ||
            categoryFilterId === CATEGORY_FILTER_ID.CELL_TYPE) &&
          !isFilterEnabled
        );
      })
      .reduce((accum, categoryFilterId: CATEGORY_FILTER_ID) => {
        accum.add(categoryFilterId);
        return accum;
      }, new Set<CATEGORY_FILTER_ID>());
  }, [isFilterEnabled]);

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

  // Store latest filter state.
  useEffect(() => {
    storeFilters(filters);
  }, [filters, storeFilters]);

  // Set row count.
  useEffect(() => {
    if (rows && rows.length) {
      setRowCount(rows.length);
    }
  }, [rows]);

  // Handle side bar open/closed state beyond scope of component.
  const [isSideBarOpen, storeIsSideBarOpen] = useSessionStorage<boolean>(
    KEYS.SIDE_BAR_COLLECTIONS,
    true
  );

  return (
    <>
      <Head>
        <title>CELL&times;GENE | Collections</title>
      </Head>
      {isError || isLoading ? null : (
        <>
          <SideBar
            label="Filters"
            isOpen={isSideBarOpen}
            onToggle={storeIsSideBarOpen}
          >
            <Filter {...filterInstance} />
          </SideBar>
          <View>
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
