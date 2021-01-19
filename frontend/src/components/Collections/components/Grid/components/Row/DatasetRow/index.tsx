import { Button, Intent, Radio } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import loadable from "@loadable/component";
import React, { FC, useEffect, useState } from "react";
import { CancelledError, QueryCache, useQueryCache } from "react-query";
import { Dataset, VALIDATION_STATUS } from "src/common/entities";
import {
  useDatasetStatus,
  useDeleteDataset,
  USE_DATASET_STATUS,
} from "src/common/queries/datasets";
import { aggregateDatasetsMetadata } from "src/components/Collections/components/Grid/common/utils";
import {
  DetailsCell,
  RightAlignedDetailsCell,
  StyledRow,
} from "src/components/Collections/components/Grid/components/Row/common/style.ts";
import { UploadingFile } from "src/components/DropboxChooser";
import DatasetUploadToast from "src/views/Collection/components/DatasetUploadToast";
import CellCount from "./components/CellCount";
import Popover from "./components/Popover";
import UploadStatus from "./components/UploadStatus";
import { TitleContainer } from "./style";
import {
  checkIfCancelled,
  checkIfFailed,
  checkIfLoading,
  FailReturn,
} from "./utils";

const FETCH_COLLECTION_INTERVAL_MS = 5 * 1000;

const AsyncTooltip = loadable(
  () =>
    /*webpackChunkName: 'Grid/Row/DatasetRow/Tooltip' */ import(
      "./components/Tooltip"
    )
);

const ErrorTooltip = ({ isFailed, error }: FailReturn) => {
  if (!isFailed) return null;

  return <AsyncTooltip error={error} />;
};

interface Props {
  dataset: Dataset;
  file?: UploadingFile;
  invalidateCollectionQuery: () => void;
  onSelect: (id: Dataset["id"]) => void;
  selected: Dataset["id"] | undefined;
}

const DatasetRow: FC<Props> = ({
  dataset,
  file,
  invalidateCollectionQuery,
  onSelect,
  selected,
}) => {
  const queryCache = useQueryCache();

  const queryResult = useDatasetStatus(
    dataset.id,
    checkIfLoading(dataset.processing_status)
  );

  const datasetStatus = queryResult.data || dataset.processing_status;

  const { upload_progress } = datasetStatus;

  const [uploadProgress, setUploadProgress] = useState(upload_progress);

  if (queryResult.isError && !(queryResult.error instanceof CancelledError)) {
    console.error(queryResult.error);
  }
  const isNamePopulated = Boolean(dataset.name);

  const name = dataset.name || file?.name || dataset.id;

  // (thuang): We need to poll the collection until the name is populated,
  // which indicates other metadata are populated too
  useCheckCollectionPopulated({
    invalidateCollectionQuery,
    isNamePopulated,
    upload_progress,
  });

  const hasFailed = checkIfFailed(datasetStatus);
  const { isFailed, error } = hasFailed;

  useCancelDatasetStatusQuery({
    datasetId: dataset.id,
    isFailed,
    isNamePopulated,
    queryCache,
  });

  useUploadProgress({
    newUploadProgress: upload_progress,
    setUploadProgress,
    uploadProgress,
  });

  const [deleteDataset] = useDeleteDataset(dataset.collection_id);

  if (checkIfCancelled(datasetStatus)) return null;

  const isLoading = checkIfLoading(datasetStatus);

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
          <Radio
            onChange={() => onSelect(dataset.id)}
            checked={selected === dataset.id}
          />
          <div>{name}</div>
          <ErrorTooltip isFailed={isFailed} error={error} />
        </TitleContainer>
        {isLoading && (
          <UploadStatus
            isValidating={
              datasetStatus.validation_status === VALIDATION_STATUS.VALIDATING
            }
            progress={datasetStatus.upload_progress}
            cancelUpload={() => {
              deleteDataset(dataset.id);
            }}
          />
        )}
      </DetailsCell>
      <Popover values={tissue} isLoading={isLoading} isFailed={isFailed} />
      <Popover values={assay} isLoading={isLoading} isFailed={isFailed} />
      <Popover values={disease} isLoading={isLoading} isFailed={isFailed} />
      <Popover values={organism} isLoading={isLoading} isFailed={isFailed} />
      <CellCount cellCount={cell_count} isLoading={isLoading} />
      <RightAlignedDetailsCell>
        {isNamePopulated && (
          <Button intent={Intent.PRIMARY} outlined text="Explore" />
        )}
      </RightAlignedDetailsCell>
    </StyledRow>
  );
};

function useCheckCollectionPopulated({
  invalidateCollectionQuery,
  isNamePopulated,
  upload_progress,
}: {
  invalidateCollectionQuery: () => void;
  isNamePopulated: boolean;
  upload_progress: number;
}) {
  useEffect(() => {
    let intervalId: number | undefined = undefined;

    if (!intervalId && upload_progress === 1 && !isNamePopulated) {
      intervalId = window?.setInterval(() => {
        if (upload_progress === 1 && !isNamePopulated) {
          invalidateCollectionQuery();
        }
      }, FETCH_COLLECTION_INTERVAL_MS);
    }

    return () => clearInterval(intervalId);
  }, [invalidateCollectionQuery, isNamePopulated, upload_progress]);
}

function useCancelDatasetStatusQuery({
  datasetId,
  isFailed,
  isNamePopulated,
  queryCache,
}: {
  datasetId: string;
  isFailed: boolean;
  isNamePopulated: boolean;
  queryCache: QueryCache;
}) {
  useEffect(() => {
    if (isFailed || isNamePopulated) {
      queryCache.cancelQueries([USE_DATASET_STATUS, datasetId]);
    }
  }, [datasetId, isFailed, isNamePopulated, queryCache]);
}

function useUploadProgress({
  newUploadProgress,
  setUploadProgress,
  uploadProgress,
}: {
  newUploadProgress: number;
  setUploadProgress: React.Dispatch<React.SetStateAction<number>>;
  uploadProgress: number;
}) {
  useEffect(() => {
    if (uploadProgress === newUploadProgress) return;

    setUploadProgress(newUploadProgress);

    if (newUploadProgress === 1) {
      DatasetUploadToast.show({
        icon: IconNames.TICK,
        intent: Intent.SUCCESS,
        message:
          "Upload was successful. Your file is being processed which will continue in the background, even if you close this window.",
      });
    }
  }, [newUploadProgress, setUploadProgress, uploadProgress]);
}

export default DatasetRow;
