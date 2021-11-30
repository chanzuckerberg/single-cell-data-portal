import Head from "next/head";
import React, { FC, useMemo } from "react";
import { CellProps, Column, useFilters, useTable } from "react-table";
import { ROUTES } from "src/common/constants/routes";
import { FilterableDataset, Ontology } from "src/common/entities";
import { FEATURES } from "src/common/featureFlags/features";
import {
  CATEGORY_KEY,
  OntologyCategoryKey,
  useFacetedFilter,
} from "src/common/hooks/useFacetedFilter";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import Categories from "src/components/Categories";
import DatasetsGrid from "src/components/Datasets/components/Grid/components/DatasetsGrid";
import filterableDatasets from "../../../tests/features/fixtures/datasets/filterable-datasets";

const Datasets: FC = () => {
  // Column configuration backing table.
  const columnConfig: Column<FilterableDataset>[] = useMemo(
    () => [
      // TODO(cc) required for "group by" of counts
      {
        Header: "Dataset ID",
        accessor: "id",
      },
      {
        Cell: DatasetNameCell,
        Header: "Dataset",
        accessor: "name", // TODO(cc) collection name
      },
      // Hidden, required for accessing collection name. TODO(cc) revisit - is this necessary?
      {
        Header: "Collection Name",
        accessor: "collection_name",
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
          "id", // TODO(cc) constant
          "collection_name", // TODO(cc) constant
          // CATEGORY_KEY.CELL_TYPE,
          // CATEGORY_KEY.IS_PRIMARY_DATA,
          // CATEGORY_KEY.SEX,
        ],
      },
    },
    useFilters
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
};

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

export default Datasets;
