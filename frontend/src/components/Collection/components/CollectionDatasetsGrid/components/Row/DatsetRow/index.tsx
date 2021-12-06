/* Copied from src/components/Collections/components/Grid/components/Row/DatasetRow and modified as an intermediate 
   upgrade to the collection datasets table. Ideally collection datasets table should be moved to react-table but
   requires changes to how dataset information (statuses etc) is calculated and rolled out. */
import { Intent, Tooltip } from "@blueprintjs/core";
import loadable from "@loadable/component";
import { FC } from "react";
import { CancelledError, useQueryCache } from "react-query";
import { PLURALIZED_METADATA_LABEL } from "src/common/constants/metadata";
import {
  ACCESS_TYPE,
  Collection,
  CONVERSION_STATUS,
  Dataset,
  UPLOAD_STATUS,
  VALIDATION_STATUS,
  VISIBILITY_TYPE,
} from "src/common/entities";
import { useDatasetStatus } from "src/common/queries/datasets";
import DownloadDataset from "src/components/Collections/components/Dataset/components/DownloadDataset";
import { aggregateDatasetsMetadata } from "src/components/Collections/components/Grid/common/utils";
import MoreDropdown from "src/components/Collections/components/Grid/components/Row/DatasetRow/components/MoreDropdown";
import RevisionStatusTag from "src/components/Collections/components/Grid/components/Row/DatasetRow/components/RevisionStatusTag";
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
} from "src/components/Collections/components/Grid/components/Row/DatasetRow/utils";
import { OVER_MAX_CELL_COUNT_TOOLTIP } from "src/components/common/Grid/common/constants";
import { checkIsOverMaxCellCount } from "src/components/common/Grid/common/utils";
import ActionButton from "src/components/common/Grid/components/ActionButton";
import ActionsCell from "src/components/common/Grid/components/ActionsCell";
import CountCell from "src/components/common/Grid/components/CountCell";
import NTagCell from "src/components/common/Grid/components/NTagCell";
import { RightAlignCell } from "src/components/common/Grid/components/RightAlignCell";
import { StatusTags } from "src/components/common/Grid/components/StatusTags";
import { DownloadButton } from "src/components/Datasets/components/Grid/common/utils";
import DatasetNameCell from "src/components/Datasets/components/Grid/components/DatasetNameCell";
import { UploadingFile } from "src/components/DropboxChooser";
import { Props as ChooserProps } from "src/components/DropboxChooser/index";
import exploreSVG from "/src/common/images/explore-blue.svg";

const AsyncTooltip = loadable(
  () =>
    /*webpackChunkName: 'Grid/Row/DatasetRow/Tooltip' */ import(
      "src/components/Collections/components/Grid/components/Row/DatasetRow/components/Tooltip"
    )
);

const AsyncUploadStatus = loadable(
  () =>
    /*webpackChunkName: 'Grid/Row/DatasetRow/UploadStatus' */ import(
      "src/components/Collections/components/Grid/components/Row/DatasetRow/components/UploadStatus"
    )
);

const ErrorTooltip = ({ isFailed, error, type }: FailReturn) => {
  if (!isFailed) return null;

  return <AsyncTooltip error={error} type={type} />;
};

interface Props {
  dataset: Dataset;
  file?: UploadingFile;
  invalidateCollectionQuery: () => void;
  visibility: Collection["visibility"];
  accessType?: Collection["access_type"];
  revisionsEnabled: boolean;
  onUploadFile: ChooserProps["onUploadFile"];
}

/* eslint-disable sonarjs/cognitive-complexity */
const DatasetRow: FC<Props> = ({
  dataset,
  file,
  invalidateCollectionQuery,
  visibility,
  accessType,
  revisionsEnabled,
  onUploadFile,
}) => {
  const queryCache = useQueryCache();

  const datasetStatusResult = useDatasetStatus(
    dataset.id,
    checkIfLoading(dataset.processing_status)
  );

  const datasetStatus = datasetStatusResult.data || dataset.processing_status;

  const initProgress = dataset?.processing_status?.upload_progress;

  const { upload_progress } = datasetStatus;

  if (
    datasetStatusResult.isError &&
    !(datasetStatusResult.error instanceof CancelledError)
  ) {
    console.error(datasetStatusResult.error);
  }

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
    queryCache,
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

  const isRDSSkipped = datasetStatus.rds_status === CONVERSION_STATUS.SKIPPED;

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
        <NTagCell label={PLURALIZED_METADATA_LABEL.DISEASE} values={disease} />
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
                collectionId={dataset.collection_id}
                datasetId={dataset.id}
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
            // isRDSSkipped is drilled 3 components down to `frontend/src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/DataFormat/index.tsx`
            isRDSSkipped={isRDSSkipped}
          />
          {hasCXGFile(dataset) && (
            <Tooltip
              boundary="viewport"
              content={
                isOverMaxCellCount ? OVER_MAX_CELL_COUNT_TOOLTIP : "Explore"
              }
              disabled={dataset.tombstone}
              intent={isOverMaxCellCount ? Intent.DANGER : undefined}
            >
              <ActionButton
                data-test-id="view-dataset-link"
                imageProps={exploreSVG}
                isDisabled={dataset.tombstone || isOverMaxCellCount}
                // @ts-expect-error -- revisit rel typing
                rel="noopener"
                target="_blank"
                url={dataset?.dataset_deployments[0]?.url}
              />
            </Tooltip>
          )}
        </ActionsCell>
      </td>
    </tr>
  );
};

export default DatasetRow;
