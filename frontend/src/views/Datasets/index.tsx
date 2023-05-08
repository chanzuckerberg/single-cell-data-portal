import Head from "next/head";
import React, { useEffect, useMemo } from "react";
import {
  CellProps,
  Column,
  Filters,
  Renderer,
  useFilters,
  useSortBy,
  useTable,
} from "react-table";
import { PLURALIZED_METADATA_LABEL } from "src/common/constants/metadata";
import { useCategoryFilter } from "src/common/hooks/useCategoryFilter/useCategoryFilter";
import { useExplainNewTab } from "src/common/hooks/useExplainNewTab";
import { useSessionStorage } from "src/common/hooks/useSessionStorage";
import { useFetchDatasetRows } from "src/common/queries/filter";
import { KEYS } from "src/common/sessionStorage/set";
import BottomBanner from "src/components/BottomBanner";
import Filter from "src/components/common/Filter";
import {
  CATEGORY_FILTER_ID,
  CellPropsValue,
  DatasetRow,
  MultiPanelSelectedUIState,
  RowPropsValue,
} from "src/components/common/Filter/common/entities";
import { ontologyLabelCellAccessorFn } from "src/components/common/Filter/common/utils";
import { buildTableCountSummary } from "src/components/common/Grid/common/utils";
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
import { ALIGNMENT } from "src/components/common/Grid/common/entities";

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
        Header: "Datasets",
        accessor: DATASET_NAME,
        showCountAndTotal: true,
      },
      {
        Cell: (({ value }: CellPropsValue<string[]>) => (
          <NTagCell label={PLURALIZED_METADATA_LABEL.TISSUE} values={value} />
        )) as Renderer<CellProps<DatasetRow>>,
        Header: "Tissue",
        accessor: ontologyLabelCellAccessorFn("tissue"),
        id: "tissue",
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
          <NTagCell label={PLURALIZED_METADATA_LABEL.ASSAY} values={value} />
        ),
        Header: "Assay",
        accessor: ontologyLabelCellAccessorFn("assay"),
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.ASSAY,
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
      {
        Cell: ({ value }: CellPropsValue<number | null>) => (
          <RightAlignCell>
            <CountCell cellCount={value || 0} />
          </RightAlignCell>
        ),
        Header: "Cells",
        accessor: "cell_count",
        alignment: ALIGNMENT.RIGHT,
        filter: "between",
        id: CATEGORY_FILTER_ID.CELL_COUNT,
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
        accessor: RECENCY,
        id: COLUMN_ID_RECENCY,
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
        accessor: "cellTypeCalculated",
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.CELL_TYPE_CALCULATED,
      },
      // Hidden, required for filter.
      {
        accessor: ontologyLabelCellAccessorFn("self_reported_ethnicity"),
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.SELF_REPORTED_ETHNICITY,
      },
      {
        accessor: "development_stage_ancestors",
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
      },
      // Hidden, required for filter.
      {
        accessor: "mean_genes_per_cell",
        filter: "between",
        id: CATEGORY_FILTER_ID.GENE_COUNT,
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
        accessor: "suspension_type",
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.SUSPENSION_TYPE,
      },
      // Hidden, required for filter.
      {
        accessor: "tissueCalculated",
        filter: "includesSome",
        id: CATEGORY_FILTER_ID.TISSUE_CALCULATED,
      },
    ],
    []
  );

  // Handle initial filter state and save of filter state beyond component scope.
  const [initialFilters, storeFilters] = useSessionStorage<Filters<DatasetRow>>(
    KEYS.FILTER_DATASETS,
    []
  );

  // Handle initial multi-panel filter UI state and save of multi-panel filter UI state beyond component scope.
  const [initialMultiPanelSelectedUIState, storeMultiPanelSelectedUIState] =
    useSessionStorage<MultiPanelSelectedUIState>(
      KEYS.FILTER_DATASETS_SELECTED_UI,
      {}
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
          CATEGORY_FILTER_ID.CELL_TYPE_CALCULATED,
          CATEGORY_FILTER_ID.SELF_REPORTED_ETHNICITY,
          CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
          CATEGORY_FILTER_ID.GENE_COUNT,
          CATEGORY_FILTER_ID.PUBLICATION_AUTHORS,
          CATEGORY_FILTER_ID.PUBLICATION_DATE_VALUES,
          CATEGORY_FILTER_ID.SEX,
          CATEGORY_FILTER_ID.SUSPENSION_TYPE,
          CATEGORY_FILTER_ID.TISSUE_CALCULATED,
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
  const categories = useMemo<Set<CATEGORY_FILTER_ID>>(() => {
    return Object.values(CATEGORY_FILTER_ID).reduce(
      (accum, categoryFilterId: CATEGORY_FILTER_ID) => {
        accum.add(categoryFilterId);
        return accum;
      },
      new Set<CATEGORY_FILTER_ID>()
    );
  }, []);

  // Set up filter instance.
  const filterInstance = useCategoryFilter(
    preFilteredRows,
    categories,
    filters,
    setFilter,
    initialMultiPanelSelectedUIState
  );

  // Store latest filter and multi-panel filter UI state.
  useEffect(() => {
    storeFilters(filters);
    storeMultiPanelSelectedUIState(filterInstance.multiPanelSelectedUIState);
  }, [
    filters,
    filterInstance.multiPanelSelectedUIState,
    storeFilters,
    storeMultiPanelSelectedUIState,
  ]);

  // Handle side bar open/closed state beyond scope of component.
  const [isSideBarOpen, storeIsSideBarOpen] = useSessionStorage<boolean>(
    KEYS.SIDE_BAR_DATASETS,
    true
  );

  return (
    <>
      <Head>
        <title>CELL&times;GENE | Datasets</title>
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
              <DatasetsGrid
                tableCountSummary={buildTableCountSummary(
                  rows,
                  preFilteredRows
                )}
                // @ts-expect-error -- revisit tableInstance typing
                tableInstance={tableInstance}
              />
            )}
          </View>
          <BottomBanner />
        </>
      )}
    </>
  );
}
