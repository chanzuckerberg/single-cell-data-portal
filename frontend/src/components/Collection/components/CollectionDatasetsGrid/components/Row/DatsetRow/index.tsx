import { Intent, Tooltip } from "@blueprintjs/core";
import loadable from "@loadable/component";
import { FC } from "react";
import { CancelledError, useQueryClient, UseQueryResult } from "react-query";
import { PLURALIZED_METADATA_LABEL } from "src/common/constants/metadata";
import {
  ACCESS_TYPE,
  Collection,
  CONVERSION_STATUS,
  Dataset,
  DatasetUploadStatus,
  UPLOAD_STATUS,
  VALIDATION_STATUS,
  VISIBILITY_TYPE,
} from "src/common/entities";
import { useDatasetStatus } from "src/common/queries/datasets";
import MoreDropdown from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/MoreDropdown";
import RevisionStatusTag from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/RevisionStatusTag";
import {
  checkIfCancelled,
  checkIfFailed,
  checkIfLoading,
  FailReturn,
  getConversionStatus,
  hasCXGFile,
  useCancelDatasetStatusQuery,
  useCheckCollectionFormatsPopulated,
  useCheckCollectionPopulated,
  useConversionProgress,
  useUploadProgress,
} from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/utils";
import DownloadDataset from "src/components/Collections/components/Dataset/components/DownloadDataset";
import { aggregateDatasetsMetadata } from "src/components/Collections/components/Grid/common/utils";
import { OVER_MAX_CELL_COUNT_TOOLTIP } from "src/components/common/Grid/common/constants";
import { checkIsOverMaxCellCount } from "src/components/common/Grid/common/utils";
import ActionsCell from "src/components/common/Grid/components/ActionsCell";
import CountCell from "src/components/common/Grid/components/CountCell";
import DiseaseCell from "src/components/common/Grid/components/DiseaseCell";
import NTagCell from "src/components/common/Grid/components/NTagCell";
import { RightAlignCell } from "src/components/common/Grid/components/RightAlignCell";
import { StatusTags } from "src/components/common/Grid/components/StatusTags";
import DatasetNameCell from "src/components/Datasets/components/Grid/components/DatasetNameCell";
import { UploadingFile } from "src/components/DropboxChooser";
import { Props as ChooserProps } from "src/components/DropboxChooser/index";
import { StyledPrimaryAnchorButton } from "src/components/common/Button/common/style";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import DownloadButton from "src/components/common/Grid/components/DownloadButton";

const AsyncTooltip = loadable(
  () =>
    /*webpackChunkName: 'Grid/Row/DatasetRow/Tooltip' */ import(
      "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/Tooltip"
    )
);

const AsyncUploadStatus = loadable(
  () =>
    /*webpackChunkName: 'Grid/Row/DatasetRow/UploadStatus' */ import(
      "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/UploadStatus"
    )
);

const ErrorTooltip = ({ isFailed, error, type }: FailReturn) => {
  if (!isFailed) return null;

  return <AsyncTooltip error={error} type={type} />;
};

const handlePossibleError = (
  datasetStatusResult: UseQueryResult<DatasetUploadStatus, unknown>
) => {
  if (
    datasetStatusResult.isError &&
    !(datasetStatusResult.error instanceof CancelledError)
  ) {
    console.error(datasetStatusResult.error);
  }
};

interface Props {
  collectionId: Collection["id"];
  dataset: Dataset;
  file?: UploadingFile;
  invalidateCollectionQuery: () => void;
  visibility: Collection["visibility"];
  accessType?: Collection["access_type"];
  revisionsEnabled: boolean;
  onUploadFile: ChooserProps["onUploadFile"];
}

const DatasetRow: FC<Props> = ({
  collectionId,
  dataset,
  file,
  invalidateCollectionQuery,
  visibility,
  accessType,
  revisionsEnabled,
  onUploadFile,
}) => {
  const queryClient = useQueryClient();

  const datasetStatusResult = useDatasetStatus(
    dataset.id,
    collectionId,
    checkIfLoading(dataset.processing_status)
  );

  const datasetStatus = datasetStatusResult.data || dataset.processing_status;

  const initProgress = dataset?.processing_status?.upload_progress;

  const { upload_progress } = datasetStatus;

  handlePossibleError(datasetStatusResult);

  const isNamePopulated = Boolean(dataset.name);

  const name = dataset.name || file?.name || dataset.id;

  // (thuang): We need to poll the collection until all the converted files
  // become available
  useCheckCollectionFormatsPopulated({
    dataset,
    datasetUploadStatus: datasetStatus,
    invalidateCollectionQuery,
  });

  const hasFailed = checkIfFailed(datasetStatus);

  const { isFailed, error, type } = hasFailed;

  // (thuang): We need to poll the collection until the name is populated,
  // which indicates other metadata are populated too
  useCheckCollectionPopulated({
    invalidateCollectionQuery,
    isFailed,
    isNamePopulated,
    validationStatus: datasetStatus.validation_status,
  });

  const isLoading = checkIfLoading(datasetStatus);

  useCancelDatasetStatusQuery({
    datasetId: dataset.id,
    isFailed,
    isLoading,
    queryClient,
  });

  useUploadProgress({
    initProgress: initProgress,
    progress: upload_progress,
  });

  useConversionProgress({
    datasetStatus: datasetStatusResult.data,
    initDatasetStatus: dataset.processing_status,
  });

  if (checkIfCancelled(datasetStatus)) return null;

  const { tissue, assay, disease, organism, cell_count } =
    aggregateDatasetsMetadata([dataset]);

  const isOverMaxCellCount = checkIsOverMaxCellCount(cell_count);

  return (
    <tr>
      <td>
        <DatasetNameCell name={name}>
          {!isLoading && (
            <ErrorTooltip isFailed={isFailed} error={error} type={type} />
          )}
          {isLoading && (
            <AsyncUploadStatus
              isWaiting={datasetStatus.upload_status === UPLOAD_STATUS.WAITING}
              isConverting={
                getConversionStatus(datasetStatus) ===
                CONVERSION_STATUS.CONVERTING
              }
              isValidating={
                datasetStatus.validation_status === VALIDATION_STATUS.VALIDATING
              }
              progress={datasetStatus.upload_progress}
            />
          )}
          <StatusTags>
            {revisionsEnabled && <RevisionStatusTag dataset={dataset} />}
          </StatusTags>
        </DatasetNameCell>
      </td>
      <td>
        <NTagCell label={PLURALIZED_METADATA_LABEL.TISSUE} values={tissue} />
      </td>
      <td>
        <DiseaseCell
          label={PLURALIZED_METADATA_LABEL.DISEASE}
          values={disease}
        />
      </td>
      <td>
        <NTagCell label={PLURALIZED_METADATA_LABEL.ASSAY} values={assay} />
      </td>
      <td>
        <NTagCell
          label={PLURALIZED_METADATA_LABEL.ORGANISM}
          values={organism}
        />
      </td>
      <td>
        <RightAlignCell>
          <CountCell cellCount={cell_count || 0} />
        </RightAlignCell>
      </td>
      <td>
        <ActionsCell>
          {visibility === VISIBILITY_TYPE.PRIVATE &&
            accessType === ACCESS_TYPE.WRITE && (
              <MoreDropdown
                collectionId={collectionId}
                datasetId={dataset.id}
                isPublished={!!dataset.published_at} // Dataset has been published.
                revisionsEnabled={revisionsEnabled}
                onUploadFile={onUploadFile}
                isLoading={isLoading}
                disabled={dataset.tombstone ?? false}
              />
            )}
          <DownloadDataset
            name={dataset?.name || ""}
            dataAssets={dataset?.dataset_assets}
            Button={DownloadButton}
            isDisabled={dataset.tombstone}
          />
          {hasCXGFile(dataset) && (
            <Tooltip
              boundary="viewport"
              content={OVER_MAX_CELL_COUNT_TOOLTIP}
              disabled={!isOverMaxCellCount}
              intent={Intent.DANGER}
            >
              <StyledPrimaryAnchorButton
                data-testid="view-dataset-link"
                disabled={dataset.tombstone || isOverMaxCellCount}
                href={dataset?.dataset_deployments[0]?.url}
                intent={Intent.PRIMARY}
                onClick={() =>
                  track(EVENTS.DATASET_EXPLORE_CLICKED, {
                    dataset_name: dataset?.name,
                  })
                }
                rel="noopener"
                target="_blank"
                text="Explore"
              />
            </Tooltip>
          )}
        </ActionsCell>
      </td>
    </tr>
  );
};

export default DatasetRow;
