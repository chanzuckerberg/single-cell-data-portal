import { AnchorButton, Intent } from "@blueprintjs/core";
import loadable from "@loadable/component";
import React, { FC } from "react";
import { CancelledError, useQueryCache } from "react-query";
import {
  CONVERSION_STATUS,
  Dataset,
  VALIDATION_STATUS,
} from "src/common/entities";
import { useDatasetStatus } from "src/common/queries/datasets";
import { aggregateDatasetsMetadata } from "src/components/Collections/components/Grid/common/utils";
import {
  DetailsCell,
  RightAlignedDetailsCell,
  StyledRow,
} from "src/components/Collections/components/Grid/components/Row/common/style.ts";
import { UploadingFile } from "src/components/DropboxChooser";
import CellCount from "./components/CellCount";
import Popover from "./components/Popover";
import { StyledRadio, TitleContainer } from "./style";
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
  onSelect: (id: Dataset["id"]) => void;
  selected?: Dataset["id"];
}

const DatasetRow: FC<Props> = ({
  dataset,
  file,
  invalidateCollectionQuery,
  onSelect,
  selected,
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

  // (thuang): We need to poll the collection until the name is populated,
  // which indicates other metadata are populated too
  useCheckCollectionPopulated({
    invalidateCollectionQuery,
    isNamePopulated,
    validationStatus: datasetStatus.validation_status,
  });

  const hasFailed = checkIfFailed(datasetStatus);
  const { isFailed, error, type } = hasFailed;

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
          <StyledRadio
            onChange={() => onSelect(dataset.id)}
            checked={selected === dataset.id}
          />
          <div>{name}</div>
          {!isLoading && (
            <ErrorTooltip isFailed={isFailed} error={error} type={type} />
          )}
        </TitleContainer>
        {isLoading && (
          <AsyncUploadStatus
            isConverting={
              getConversionStatus(datasetStatus) ===
              CONVERSION_STATUS.CONVERTING
            }
            isValidating={
              datasetStatus.validation_status === VALIDATION_STATUS.VALIDATING
            }
            progress={datasetStatus.upload_progress}
            datasetId={dataset.id}
            collectionId={dataset.collection_id}
          />
        )}
      </DetailsCell>
      <Popover values={tissue} isLoading={isMetadataLoading} />
      <Popover values={assay} isLoading={isMetadataLoading} />
      <Popover values={disease} isLoading={isMetadataLoading} />
      <Popover values={organism} isLoading={isMetadataLoading} />
      <CellCount cellCount={cell_count} isLoading={isMetadataLoading} />
      <RightAlignedDetailsCell>
        {hasCXGFile(dataset) && (
          <AnchorButton
            intent={Intent.PRIMARY}
            outlined
            text="Explore"
            href={dataset?.dataset_deployments[0]?.url}
            target="_blank"
            rel="noopener"
            data-test-id="view-dataset-link"
          />
        )}
      </RightAlignedDetailsCell>
    </StyledRow>
  );
};

export default DatasetRow;
