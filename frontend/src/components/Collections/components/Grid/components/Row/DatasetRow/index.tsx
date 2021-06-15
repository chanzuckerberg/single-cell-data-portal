import { AnchorButton, Classes, Intent, Tooltip } from "@blueprintjs/core";
import loadable from "@loadable/component";
import React, { FC } from "react";
import { CancelledError, useQueryCache } from "react-query";
import {
  ACCESS_TYPE,
  Collection,
  CONVERSION_STATUS,
  Dataset,
  UPLOAD_STATUS,
  VALIDATION_STATUS,
  VISIBILITY_TYPE,
} from "src/common/entities";
import ExploreSVG from "src/common/images/explore.svg";
import { useDatasetStatus } from "src/common/queries/datasets";
import DownloadDataset from "src/components/Collections/components/Dataset/components/DownloadDataset";
import { aggregateDatasetsMetadata } from "src/components/Collections/components/Grid/common/utils";
import {
  ActionButton,
  ActionButtonsContainer,
  ActionCell,
  DetailsCell,
  StyledRow,
} from "src/components/Collections/components/Grid/components/Row/common/style";
import { UploadingFile } from "src/components/DropboxChooser";
import { Props as ChooserProps } from "src/components/DropboxChooser/index";
import CellCount from "./components/CellCount";
import DownloadButton from "./components/DownloadButton";
import MoreDropdown from "./components/MoreDropdown";
import Popover from "./components/Popover";
import RevisionStatusTag from "./components/RevisionStatusTag";
import { TitleContainer } from "./style";
import {
  checkIfCancelled,
  checkIfFailed,
  checkIfLoading,
  checkIfMetadataLoading,
  FailReturn,
  getConversionStatus,
  hasCXGFile,
  useCancelDatasetStatusQuery,
  useCheckCollectionFormatsPopulated,
  useCheckCollectionPopulated,
  useConversionProgress,
  useUploadProgress,
} from "./utils";

const AsyncTooltip = loadable(
  () =>
    /*webpackChunkName: 'Grid/Row/DatasetRow/Tooltip' */ import(
      "./components/Tooltip"
    )
);

const AsyncUploadStatus = loadable(
  () =>
    /*webpackChunkName: 'Grid/Row/DatasetRow/UploadStatus' */ import(
      "./components/UploadStatus"
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

  const isMetadataLoading = checkIfMetadataLoading(dataset, datasetStatus);

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

  const {
    tissue,
    assay,
    disease,
    organism,
    cell_count,
  } = aggregateDatasetsMetadata([dataset]);

  return (
    <StyledRow>
      <DetailsCell>
        <TitleContainer>
          <span>{name}</span>
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
          {revisionsEnabled && <RevisionStatusTag dataset={dataset} />}
        </TitleContainer>
      </DetailsCell>
      <Popover values={tissue} isLoading={isMetadataLoading} />
      <Popover values={assay} isLoading={isMetadataLoading} />
      <Popover values={disease} isLoading={isMetadataLoading} />
      <Popover values={organism} isLoading={isMetadataLoading} />
      <CellCount cellCount={cell_count} isLoading={isMetadataLoading} />
      <ActionCell>
        <ActionButtonsContainer>
          <ActionButton>
            {visibility === VISIBILITY_TYPE.PRIVATE &&
              accessType === ACCESS_TYPE.WRITE && (
                <MoreDropdown
                  collectionId={dataset.collection_id}
                  datasetId={dataset.id}
                  revisionsEnabled={revisionsEnabled}
                  onUploadFile={onUploadFile}
                  isLoading={isLoading}
                />
              )}
          </ActionButton>

          <ActionButton>
            <DownloadDataset
              name={dataset?.name || ""}
              dataAssets={dataset?.dataset_assets}
              Button={DownloadButton}
            />
          </ActionButton>

          <ActionButton>
            {hasCXGFile(dataset) && (
              <Tooltip content="Explore">
                <AnchorButton
                  minimal
                  intent={Intent.PRIMARY}
                  icon={
                    <span className={Classes.ICON}>
                      <img alt="Explore" src={String(ExploreSVG)} />
                    </span>
                  }
                  href={dataset?.dataset_deployments[0]?.url}
                  target="_blank"
                  rel="noopener"
                  data-test-id="view-dataset-link"
                />
              </Tooltip>
            )}
          </ActionButton>
        </ActionButtonsContainer>
      </ActionCell>
    </StyledRow>
  );
};

export default DatasetRow;
