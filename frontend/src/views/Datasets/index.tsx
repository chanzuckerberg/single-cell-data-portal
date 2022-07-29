import Head from "next/head";
import React, { useEffect, useMemo } from "react";
import { Column, Filters, useFilters, useSortBy, useTable } from "react-table";
import { PLURALIZED_METADATA_LABEL } from "src/common/constants/metadata";
import { FEATURES } from "src/common/featureFlags/features";
import {
  FilterCategoryKey,
  useCategoryFilter,
} from "src/common/hooks/useCategoryFilter";
import { useExplainNewTab } from "src/common/hooks/useExplainNewTab";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { useSessionStorage } from "src/common/hooks/useSessionStorage";
import { useFetchDatasetRows } from "src/common/queries/filter";
import { KEYS } from "src/common/sessionStorage/set";
import Filter from "src/components/common/Filter";
import {
  CellPropsValue,
  DatasetRow,
  FILTER_CATEGORY_KEY,
  RowPropsValue,
} from "src/components/common/Filter/common/entities";
import { ontologyCellAccessorFn } from "src/components/common/Filter/common/utils";
import CountCell from "src/components/common/Grid/components/CountCell";
import DiseaseCell from "src/components/common/Grid/components/DiseaseCell";
import { GridHero } from "src/components/common/Grid/components/Hero";
import NTagCell from "src/components/common/Grid/components/NTagCell";
import { RightAlignCell } from "src/components/common/Grid/components/RightAlignCell";
import SideBar from "src/components/common/SideBar";
import DatasetsActionsCell from "src/components/Datasets/components/Grid/components/DatasetActionsCell";
import DatasetNameCell from "src/components/Datasets/components/Grid/components/DatasetNameCell";
import { DatasetsGrid } from "src/components/Datasets/components/Grid/components/DatasetsGrid/style";
import { View } from "../globalStyle";

/**
 * Collection ID object key.
 */
const COLLECTION_ID = "collection_id";

/**
 * Collection name object key.
 */
const COLLECTION_NAME = "collection_name";

/**
 * Dataset ID object key.
 */
const DATASET_ID = "id";

/**
 * Dataset name object key.
 */
const DATASET_NAME = "name";

/**
 * Explorer URL object key.
 */
const EXPLORER_URL = "explorer_url";

/**
 * isOverMaxCellCount object key.
 */
const IS_OVER_MAX_CELL_COUNT = "isOverMaxCellCount";

/**
 * Recency object key.
 */
const RECENCY = "recency";

/**
 * Key identifying recency sort by column.
 */
const COLUMN_ID_RECENCY = "recency";

export default function Datasets(): JSX.Element {
  /* Pop toast if user has come from Explorer with work in progress */
  useExplainNewTab(
    "To maintain your in-progress work on the previous dataset, we opened a new tab."
  );

  // Filterable datasets joined from datasets and collections responses.
  const { isError, isLoading, rows: datasetRows } = useFetchDatasetRows();

  // Column configuration backing table.
  const columnConfig: Column<DatasetRow>[] = useMemo(
    () => [
      {
        Cell: ({ row: { values } }: RowPropsValue<DatasetRow>) => {
          return (
            <DatasetNameCell
              collectionId={values.collection_id}
              collectionName={values.collection_name}
              name={values.name}
            />
          );
        },
        Header: "Dataset",
        accessor: DATASET_NAME,
      },
      {
        Cell: ({ value }: CellPropsValue<string[]>) => (
          <NTagCell label={PLURALIZED_METADATA_LABEL.TISSUE} values={value} />
        ),
        Header: "Tissue",
        accessor: ontologyCellAccessorFn("tissue"),
        filter: "includesSome",
        id: FILTER_CATEGORY_KEY.TISSUE_DEPRECATED,
      },
      {
        Cell: ({ value }: CellPropsValue<string[]>) => (
          <DiseaseCell
            label={PLURALIZED_METADATA_LABEL.DISEASE}
            values={value}
          />
        ),
        Header: "Disease",
        accessor: ontologyCellAccessorFn("disease"),
        filter: "includesSome",
        id: FILTER_CATEGORY_KEY.DISEASE,
      },
      {
        Cell: ({ value }: CellPropsValue<string[]>) => (
          <NTagCell label={PLURALIZED_METADATA_LABEL.ASSAY} values={value} />
        ),
        Header: "Assay",
        accessor: ontologyCellAccessorFn("assay"),
        filter: "includesSome",
        id: FILTER_CATEGORY_KEY.ASSAY,
      },
      {
        Cell: ({ value }: CellPropsValue<string[]>) => (
          <NTagCell label={PLURALIZED_METADATA_LABEL.ORGANISM} values={value} />
        ),
        Header: "Organism",
        accessor: ontologyCellAccessorFn("organism"),
        filter: "includesSome",
        id: FILTER_CATEGORY_KEY.ORGANISM,
      },
      {
        Cell: ({ value }: CellPropsValue<number>) => (
          <RightAlignCell>
            <CountCell cellCount={value || 0} />
          </RightAlignCell>
        ),
        Header: <RightAlignCell>Cells</RightAlignCell>,
        accessor: "cell_count",
        filter: "between",
        id: FILTER_CATEGORY_KEY.CELL_COUNT,
      },
      {
        Cell: ({ row: { values } }: RowPropsValue<DatasetRow>) => (
          <DatasetsActionsCell
            datasetId={values.id}
            isOverMaxCellCount={values.isOverMaxCellCount}
            name={values.name}
            tombstone={false} // Only public datasets are displayed in the datasets index.
            explorerUrl={values.explorer_url}
          />
        ),
        accessor: (datasetRow: DatasetRow): DatasetRow => datasetRow,
        id: "dataset_actions",
      },
      // Hidden, required for sorting by recency.
      {
        accessor: RECENCY, // TODO(cc) review need for constant here as this is restricted by type (here and collections)
        id: COLUMN_ID_RECENCY, // TODO(cc) do we need this?
      },
      // Hidden, required for accessing dataset ID via row.values, for download functionality.
      {
        accessor: DATASET_ID,
      },
      // Hidden, required for accessing collection ID via row.values, for building link to collection detail page.
      {
        accessor: COLLECTION_ID,
      },
      // Hidden, required for accessing collection name via row.values, for display.
      {
        accessor: COLLECTION_NAME,
      },
      // Hidden, required for accessing explorer_url via row.values, for display.
      {
        accessor: EXPLORER_URL,
      },
      // Hidden, required for accessing isOverMaxCellCount via row.values, for Explore functionality.
      {
        accessor: IS_OVER_MAX_CELL_COUNT,
      },
      // Hidden, required for filter.
      {
        accessor: ontologyCellAccessorFn("cell_type"),
        filter: "includesSome",
        id: FILTER_CATEGORY_KEY.CELL_TYPE_DEPRECATED,
      },
      // Hidden, required for filter.
      {
        accessor: "cell_type_ancestors",
        filter: "includesSome",
        id: FILTER_CATEGORY_KEY.CELL_TYPE,
      },
      // Hidden, required for filter.
      {
        accessor: ontologyCellAccessorFn("ethnicity"),
        filter: "includesSome",
        id: FILTER_CATEGORY_KEY.ETHNICITY,
      },
      {
        accessor: "development_stage_ancestors",
        filter: "includesSome",
        id: FILTER_CATEGORY_KEY.DEVELOPMENT_STAGE,
      },
      // Hidden, required for filter.
      {
        accessor: "mean_genes_per_cell",
        filter: "between",
        id: FILTER_CATEGORY_KEY.GENE_COUNT,
      },
      // Hidden, required for filter.
      {
        accessor: "publicationAuthors",
        filter: "includesSome",
        id: FILTER_CATEGORY_KEY.PUBLICATION_AUTHORS,
      },
      // Hidden, required for filter.
      {
        accessor: "publicationDateValues",
        filter: "includesSome",
        id: FILTER_CATEGORY_KEY.PUBLICATION_DATE_VALUES,
      },
      // Hidden, required for filter.
      {
        accessor: ontologyCellAccessorFn("sex"),
        filter: "includesSome",
        id: FILTER_CATEGORY_KEY.SEX,
      },
      // Hidden, required for filter.
      {
        accessor: "tissue_ancestors",
        filter: "includesSome",
        id: FILTER_CATEGORY_KEY.TISSUE,
      },
    ],
    []
  );

  // Handle initial filter state and save of filter state beyond component scope.
  const [initialFilters, storeFilters] = useSessionStorage<Filters<DatasetRow>>(
    KEYS.FILTER_DATASETS,
    []
  );

  // Table init
  const tableInstance = useTable<DatasetRow>(
    {
      columns: columnConfig,
      data: datasetRows,
      initialState: {
        filters: initialFilters,
        hiddenColumns: [
          DATASET_ID,
          COLLECTION_ID,
          COLLECTION_NAME,
          COLUMN_ID_RECENCY,
          FILTER_CATEGORY_KEY.CELL_TYPE_DEPRECATED,
          FILTER_CATEGORY_KEY.CELL_TYPE,
          FILTER_CATEGORY_KEY.ETHNICITY,
          FILTER_CATEGORY_KEY.DEVELOPMENT_STAGE,
          FILTER_CATEGORY_KEY.GENE_COUNT,
          FILTER_CATEGORY_KEY.PUBLICATION_AUTHORS,
          FILTER_CATEGORY_KEY.PUBLICATION_DATE_VALUES,
          FILTER_CATEGORY_KEY.SEX,
          FILTER_CATEGORY_KEY.TISSUE,
          EXPLORER_URL,
          IS_OVER_MAX_CELL_COUNT,
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

  // Determine the set of categories to display for the datasets view.
  const isFilterEnabled = useFeatureFlag(FEATURES.FILTER); // TODO(cc) remove with #2569.
  const categories = useMemo<Set<FILTER_CATEGORY_KEY>>(() => {
    return Object.values(FILTER_CATEGORY_KEY)
      .filter(
        (categoryKey: FilterCategoryKey) =>
          !(
            (categoryKey === FILTER_CATEGORY_KEY.TISSUE ||
              categoryKey === FILTER_CATEGORY_KEY.CELL_TYPE) &&
            !isFilterEnabled
          )
      )
      .reduce((accum, categoryKey: FilterCategoryKey) => {
        accum.add(categoryKey);
        return accum;
      }, new Set<FILTER_CATEGORY_KEY>());
  }, [isFilterEnabled]);

  // Set up filter instance.
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

  // Handle side bar open/closed state beyond scope of component.
  const [isSideBarOpen, storeIsSideBarOpen] = useSessionStorage<boolean>(
    KEYS.SIDE_BAR_DATASETS,
    true
  );

  return (
    <>
      <Head>
        <title>cellxgene | Datasets</title>
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
                <p>There are no datasets matching those filters.</p>
              </GridHero>
            ) : (
              // @ts-expect-error -- revisit tableInstance typing
              <DatasetsGrid tableInstance={tableInstance} />
            )}
          </View>
        </>
      )}
    </>
  );
}
