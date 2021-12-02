import Head from "next/head";
import React, { useMemo, useState } from "react";
import {
  CellProps,
  Column,
  useFilters,
  useSortBy,
  useTable,
} from "react-table";
import { ROUTES } from "src/common/constants/routes";
import { Ontology } from "src/common/entities";
import { FEATURES } from "src/common/featureFlags/features";
import {
  OntologyCategoryKey,
  useFacetedFilter,
} from "src/common/hooks/useFacetedFilter";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { fetchFilterableDatasets } from "src/common/queries/filterable-datasets";
import Categories from "src/components/Categories";
import {
  CATEGORY_KEY,
  FilterableDataset,
} from "src/components/common/Filter/common/entities";
import DatasetsGrid from "src/components/Datasets/components/Grid/components/DatasetsGrid";

// Collection name object key
const COLLECTION_NAME = "collection_name";

// Dataset ID object key
const DATASET_ID = "id";

// Dataset name object key
const DATASET_NAME = "name";

// Key identifying recency sort by column
const COLUMN_ID_RECENCY = "recency";

export default function Datasets(): JSX.Element {
  // Filterable datasets joined from datasets index and collections index responses.
  const [filterableDatasets] = useState<FilterableDataset[]>(
    fetchFilterableDatasets()
  );

  // Column configuration backing table.
  const columnConfig: Column<FilterableDataset>[] = useMemo(
    () => [
      // Hidden, ID column, required for "group by" functionality during category summarization. TODO(cc) revisit once collection functionality is confirmed.
      {
        accessor: DATASET_ID,
      },
      // Hidden, required for sorting by recency.
      {
        accessor: (dataset: FilterableDataset): number =>
          dataset.revised_at ?? dataset.published_at,
        id: COLUMN_ID_RECENCY,
      },
      // Hidden, required for accessing collection name via row.values, for display.
      {
        accessor: COLLECTION_NAME,
      },
      {
        Cell: DatasetNameCell,
        Header: "Dataset",
        accessor: DATASET_NAME,
      },
      {
        Cell: Cell,
        Header: "Tissue",
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.TISSUE),
        filter: "includesSome",
        id: CATEGORY_KEY.TISSUE,
      },
      {
        Cell: Cell,
        Header: "Disease",
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.DISEASE),
        filter: "includesSome",
        id: CATEGORY_KEY.DISEASE,
      },
      {
        Cell: Cell,
        Header: "Assay",
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.ASSAY),
        filter: "includesSome",
        id: CATEGORY_KEY.ASSAY,
      },
      {
        Cell: Cell,
        Header: "Organism",
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.ORGANISM),
        filter: "includesSome",
        id: CATEGORY_KEY.ORGANISM,
      },
      {
        Header: "Cells",
        accessor: "cell_count",
      },
      {
        Cell: Cell,
        Header: "Cell Type",
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.CELL_TYPE),
        filter: "includesSome",
        id: CATEGORY_KEY.CELL_TYPE,
      },
      {
        Header: "Primary Data",
        accessor: CATEGORY_KEY.IS_PRIMARY_DATA,
        filter: "includesSome",
      },
      {
        Cell: Cell,
        Header: "Sex",
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.SEX),
        filter: "includesSome",
        id: CATEGORY_KEY.SEX,
      },
    ],
    []
  );

  // Table init
  const tableInstance = useTable<FilterableDataset>(
    {
      columns: columnConfig,
      data: filterableDatasets,
      initialState: {
        hiddenColumns: [
          DATASET_ID,
          COLLECTION_NAME,
          COLUMN_ID_RECENCY,
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
    "id"
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
        <title>cellxgene | Datasets</title>
      </Head>
      <div style={{ display: "flex" }}>
        <Categories {...filterInstance} />
        <DatasetsGrid tableInstance={tableInstance} />
      </div>
    </>
  );
}

/**
 * Generic table cell component displaying multi-value cell values.
 * @param props - Cell-specific properties supplied from react-table.
 * @returns Array of DOM elements, one for each value in multi-value cell.
 */
function Cell(props: CellProps<FilterableDataset, string[]>): JSX.Element[] {
  const {
    cell: { value },
  } = props;
  return value.map((v: string) => <div key={v}>{v}</div>);
}

/**
 * Table cell component displaying dataset and collection names.
 * @param props - Cell-specific properties supplied from react-table.
 * @returns DOM element containing both the dataset and collection names.
 */
function DatasetNameCell(
  props: CellProps<FilterableDataset, string[]>
): JSX.Element {
  return (
    <div>
      <div>{props.row.values.name}</div>
      <div style={{ color: "#5C7080", fontSize: "12px" }}>
        {props.row.values.collection_name}
      </div>
    </div>
  );
}

/**
 * Create function to be used by column.accessor in react-table column definition, for columns containing ontology
 * metadata (ontology label and key) values.
 * @param key - Object key of value to display in cell.
 * @returns Function that returns value with the given key, to display in a cell.
 */
function ontologyCellAccessorFn(key: OntologyCategoryKey) {
  return (dataset: FilterableDataset) =>
    dataset[key].map((o: Ontology) => o.label);
}
