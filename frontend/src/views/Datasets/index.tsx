import { IButtonProps, Intent, Tooltip } from "@blueprintjs/core";
import Head from "next/head";
import React, { useMemo } from "react";
import {
  CellProps,
  Column,
  useFilters,
  useSortBy,
  useTable,
} from "react-table";
import { PLURALIZED_METADATA_LABEL } from "src/common/constants/metadata";
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
} from "src/components/common/Filter/common/entities";
import { ontologyCellAccessorFn } from "src/components/common/Filter/common/utils";
import ActionButton from "src/components/common/Grid/components/ActionButton";
import ActionsCell from "src/components/common/Grid/components/ActionsCell";
import CountCell from "src/components/common/Grid/components/CountCell";
import { GridHero } from "src/components/common/Grid/components/Hero";
import LinkCell from "src/components/common/Grid/components/LinkCell";
import NTagCell from "src/components/common/Grid/components/NTagCell";
import { RightAlignCell } from "src/components/common/Grid/components/RightAlignCell";
import { SubTitle } from "src/components/common/Grid/components/SubTitle";
import { Title } from "src/components/common/Grid/components/Title";
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

// TODO(cc) temp hard coded - remove once dataset_assets is returned from /datasets/index
const DATASET_ASSETS = [
  {
    // eslint-disable-next-line sonarjs/no-duplicate-string -- will be removed as per above
    dataset_id: "aced2b9a-0107-4b06-9dae-059d170b94a6",
    filename: "local.h5ad",
    filetype: "H5AD",
    id: "fcb3a516-d66d-4825-8dac-48a3b7788404",
    s3_uri:
      "s3://corpora-data-staging/14826f8b-2bdf-492d-bc35-928580b83d94/local.h5ad",
    type: "REMIX",
  },
  {
    dataset_id: "aced2b9a-0107-4b06-9dae-059d170b94a6",
    filename: "local.rds",
    filetype: "RDS",
    id: "7595110b-dee6-40ed-bf46-b4365793c2a0",
    s3_uri:
      "s3://corpora-data-staging/14826f8b-2bdf-492d-bc35-928580b83d94/local.rds",
    type: "REMIX",
  },
  {
    dataset_id: "aced2b9a-0107-4b06-9dae-059d170b94a6",
    filename: "explorer_cxg",
    filetype: "CXG",
    id: "9529083a-6bf6-40d8-aaf0-d1e6ba541727",
    s3_uri:
      "s3://hosted-cellxgene-staging/14826f8b-2bdf-492d-bc35-928580b83d94.cxg/",
    type: "REMIX",
  },
] as DatasetAsset[];

// TODO(cc) temp hard coded - remove once dataset_assets is returned from /datasets/index
const DATASET_DEPLOYMENTS = [
  {
    url: "https://cellxgene.staging.single-cell.czi.technology/e/aced2b9a-0107-4b06-9dae-059d170b94a6.cxg/",
  },
];

export default function Datasets(): JSX.Element {
  // Filterable datasets joined from datasets and collections responses.
  const { isError, isLoading, rows: datasetRows } = useFetchDatasetRows();

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
      {isError || isLoading ? null : (
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
              // @ts-expect-error -- revisit tableInstance typing
              <DatasetsGrid tableInstance={tableInstance} />
            )}
          </View>
        </>
      )}
    </>
  );
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
  const DownloadButton = (downloadProps: IButtonProps) => (
    <ActionButton imageProps={downloadSVG} {...downloadProps} />
  );
  const isOverMaxCellCount = checkIsOverMaxCellCount(cell_count);

  return (
    <ActionsCell>
      <DownloadDataset
        Button={DownloadButton}
        dataAssets={DATASET_ASSETS} // TODO(cc) dataset_assets
        isDisabled={false} // tombstone will always be false
        isRDSSkipped={false} // TODO(cc)
        name={name}
      />
      <Tooltip // hasCXGFile(dataset) TODO(cc)
        content={isOverMaxCellCount ? OVER_MAX_CELL_COUNT_TOOLTIP : "Explore"}
        intent={isOverMaxCellCount ? Intent.DANGER : undefined}
      >
        <ActionButton
          data-test-id="view-dataset-link" // TODO(cc)
          imageProps={exploreSVG}
          isDisabled={isOverMaxCellCount}
          // @ts-expect-error -- revisit rel typing
          rel="noopener"
          target="_blank"
          url={DATASET_DEPLOYMENTS[0].url} // dataset?.dataset_deployments[0]?.url
        />
      </Tooltip>
    </ActionsCell>
  );
}

/**
 * Table cell component displaying dataset and collection names.
 * @param props - Cell-specific properties supplied from react-table.
 * @returns DOM element containing both the dataset and collection names.
 */
function DatasetNameCell(props: CellProps<DatasetRow, string[]>): JSX.Element {
  const {
    row: { values },
  } = props;
  const { collection_id, collection_name, name } = values;
  const url = ROUTES.COLLECTION.replace(":id", collection_id);
  return (
    <>
      <Title>{name}</Title>
      {collection_id ? (
        <SubTitle>
          <LinkCell url={url} value={collection_name} />
        </SubTitle>
      ) : (
        <SubTitle>{collection_name}</SubTitle>
      )}
    </>
  );
}
