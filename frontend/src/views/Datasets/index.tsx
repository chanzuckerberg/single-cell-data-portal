import { Intent, Tooltip } from "@blueprintjs/core";
import Head from "next/head";
import React, { useMemo } from "react";
import {
  CellProps,
  Column,
  useFilters,
  useSortBy,
  useTable,
} from "react-table";
import { ROUTES } from "src/common/constants/routes";
import { DatasetAsset } from "src/common/entities";
import { FEATURES } from "src/common/featureFlags/features";
import { useCategoryFilter } from "src/common/hooks/useCategoryFilter";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { useFetchDatasetRows } from "src/common/queries/filter";
import DownloadDataset from "src/components/Collections/components/Dataset/components/DownloadDataset";
import {
  checkIsOverMaxCellCount,
  OVER_MAX_CELL_COUNT_TOOLTIP,
} from "src/components/Collections/components/Grid/components/Row/DatasetRow";
import Filter from "src/components/common/Filter";
import {
  CATEGORY_KEY,
  CellPropsValue,
  DatasetRow,
  PLURALIZED_METADATA_LABEL,
} from "src/components/common/Filter/common/entities";
import { ontologyCellAccessorFn } from "src/components/common/Filter/common/utils";
import CellActionButton from "src/components/common/Grid/components/CellActionButton";
import CellActions from "src/components/common/Grid/components/CellActions";
import CountCell from "src/components/common/Grid/components/CountCell";
import { GridHero } from "src/components/common/Grid/components/Hero";
import NTagCell from "src/components/common/Grid/components/NTagCell";
import { RightAlignCell } from "src/components/common/Grid/components/RightAlignCell";
import SideBar from "src/components/common/SideBar";
import { DatasetsGrid } from "src/components/Datasets/components/Grid/components/DatasetsGrid/style";
import { View } from "../globalStyle";
import downloadSVG from "/src/common/images/download-blue.svg";
import exploreSVG from "/src/common/images/explore-blue.svg";

// Collection ID object key.
const COLLECTION_ID = "collection_id";

// Collection name object key.
const COLLECTION_NAME = "collection_name";

// Dataset ID object key.
const DATASET_ID = "id";

// Dataset name object key.
const DATASET_NAME = "name";

// Key identifying recency sort by column.
const COLUMN_ID_RECENCY = "recency";

const DATASET_ASSETS = [
  {
    created_at: 1637192748.443175,
    dataset_id: "218acb0f-9f2f-4f76-b90b-15a4b7c7f629",
    filename: "explorer_cxg",
    filetype: "CXG",
    id: "d5b661e5-c71b-45f9-8ac0-b54b41b263e1",
    s3_uri:
      "s3://hosted-cellxgene-prod/218acb0f-9f2f-4f76-b90b-15a4b7c7f629.cxg/",
    type: "REMIX",
    updated_at: 1637192748.443183,
    user_submitted: true,
  },
  {
    created_at: 1637192883.138201,
    dataset_id: "218acb0f-9f2f-4f76-b90b-15a4b7c7f629",
    filename: "local.h5ad",
    filetype: "H5AD",
    id: "b51bfa2d-22b2-4c65-9803-d36d4de973fa",
    s3_uri:
      "s3://corpora-data-prod/218acb0f-9f2f-4f76-b90b-15a4b7c7f629/local.h5ad",
    type: "REMIX",
    updated_at: 1637192883.138206,
    user_submitted: true,
  },
] as unknown as DatasetAsset[];

export default function Datasets(): JSX.Element {
  // Filterable datasets joined from datasets and collections responses.
  const { error, loading, rows: datasetRows } = useFetchDatasetRows();

  // Column configuration backing table.
  const columnConfig: Column<DatasetRow>[] = useMemo(
    () => [
      {
        Cell: DatasetNameCell,
        Header: "Dataset",
        accessor: DATASET_NAME,
      },
      {
        Cell: ({ value }: CellPropsValue) => (
          <NTagCell label={PLURALIZED_METADATA_LABEL.TISSUE} values={value} />
        ),
        Header: "Tissue",
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.TISSUE),
        filter: "includesSome",
        id: CATEGORY_KEY.TISSUE,
      },
      {
        Cell: ({ value }: CellPropsValue) => (
          <NTagCell label={PLURALIZED_METADATA_LABEL.DISEASE} values={value} />
        ),
        Header: "Disease",
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.DISEASE),
        filter: "includesSome",
        id: CATEGORY_KEY.DISEASE,
      },
      {
        Cell: ({ value }: CellPropsValue) => (
          <NTagCell label={PLURALIZED_METADATA_LABEL.ASSAY} values={value} />
        ),
        Header: "Assay",
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.ASSAY),
        filter: "includesSome",
        id: CATEGORY_KEY.ASSAY,
      },
      {
        Cell: ({ value }: CellPropsValue) => (
          <NTagCell label={PLURALIZED_METADATA_LABEL.ORGANISM} values={value} />
        ),
        Header: "Organism",
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.ORGANISM),
        filter: "includesSome",
        id: CATEGORY_KEY.ORGANISM,
      },
      {
        Cell: ({ value }: CellPropsValue) => (
          <RightAlignCell>
            <CountCell value={value} />
          </RightAlignCell>
        ),
        Header: <RightAlignCell>Cells</RightAlignCell>,
        accessor: "cell_count",
      },
      {
        Cell: DatasetsActionsCell,
        accessor: (datasetRow: DatasetRow): DatasetRow => datasetRow,
        id: "dataset_actions",
      },
      // Hidden, required for sorting by recency.
      {
        accessor: (datasetRow: DatasetRow): number =>
          datasetRow.revised_at ?? datasetRow.published_at,
        id: COLUMN_ID_RECENCY,
      },
      // Hidden, required for accessing collection ID via row.values, for building link to collection detail page.
      {
        accessor: COLLECTION_ID,
      },
      // Hidden, required for accessing collection name via row.values, for display.
      {
        accessor: COLLECTION_NAME,
      },
      // Hidden, required for filter.
      {
        accessor: ontologyCellAccessorFn(CATEGORY_KEY.CELL_TYPE),
        filter: "includesSome",
        id: CATEGORY_KEY.CELL_TYPE,
      },
      // Hidden, required for filter.
      {
        accessor: CATEGORY_KEY.IS_PRIMARY_DATA,
        filter: "includesSome",
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
  const tableInstance = useTable<DatasetRow>(
    {
      columns: columnConfig,
      data: datasetRows,
      initialState: {
        hiddenColumns: [
          DATASET_ID,
          COLLECTION_ID,
          COLLECTION_NAME,
          COLUMN_ID_RECENCY,
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
        <title>cellxgene | Datasets</title>
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
                <p>There are no datasets matching those filters.</p>
              </GridHero>
            ) : (
              <DatasetsGrid tableInstance={tableInstance} />
            )}
          </View>
        </>
      )}
    </>
  );
}

/**
 * Returns action button with corresponding action icon.
 * @param imageData
 * @returns cell action button with corresponding action icon
 */
function ActionButton(imageData: StaticImageData): JSX.Element {
  return <CellActionButton imageProps={imageData} />;
}

/**
 * Table cell component displaying dataset actions "download" and "explore".
 * @param props - Cell-specific properties supplied from react-table.
 * @returns DOM element containing dataset actions "download" and "explore".
 */
function DatasetsActionsCell(
  props: CellProps<DatasetRow, string[]>
): JSX.Element {
  const {
    row: { values },
  } = props;
  const { cell_count, name } = values || "";
  const DownloadButton = () => ActionButton(downloadSVG);
  const isOverMaxCellCount = checkIsOverMaxCellCount(cell_count);
  // hasCXGFile(dataset) TODO(cc) test

  return (
    <CellActions>
      <DownloadDataset
        Button={DownloadButton}
        dataAssets={DATASET_ASSETS} // TODO(cc) dataset_assets
        isDisabled={false} // tombstone will always be false
        isRDSSkipped={false} // TODO(cc)
        name={name}
      />
      <Tooltip
        content={isOverMaxCellCount ? OVER_MAX_CELL_COUNT_TOOLTIP : "Explore"}
        intent={isOverMaxCellCount ? Intent.DANGER : undefined}
      >
        <CellActionButton
          data-test-id="view-dataset-link" // TODO(cc)
          imageProps={exploreSVG}
          isDisabled={true} // dataset.tombstone || isOverMaxCellCount
          rel="noopener"
          target="_blank"
          url={"/"} // dataset?.dataset_deployments[0]?.url
        />
      </Tooltip>
    </CellActions>
  );
}

/**
 * Table cell component displaying dataset and collection names.
 * @param props - Cell-specific properties supplied from react-table.
 * @returns DOM element containing both the dataset and collection names.
 */
function DatasetNameCell(props: CellProps<DatasetRow, string[]>): JSX.Element {
  return (
    <div>
      <div>{props.row.values.name}</div>
      <div style={{ color: "#5C7080", fontSize: "12px" }}>
        {props.row.values.collection_name}
      </div>
    </div>
  );
}
