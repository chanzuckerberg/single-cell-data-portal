import { Intent } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { useEffect } from "react";
import { QueryClient } from "react-query";
import {
  CONVERSION_STATUS,
  Dataset,
  DatasetUploadStatus,
  PROCESSING_STATUS,
  UPLOAD_STATUS,
  VALIDATION_STATUS,
} from "src/common/entities";
import { USE_DATASET_STATUS } from "src/common/queries/datasets";
import Toast from "src/views/Collection/components/Toast";

export enum FAILED_RETURN_TYPE {
  UPLOAD = "UPLOAD",
  VALIDATION = "VALIDATION",
  CONVERSION = "CONVERSIOn",
  PROCESS = "PROCESS",
}

export type FailReturn = {
  isFailed: boolean;
  error?:
    | VALIDATION_STATUS
    | UPLOAD_STATUS
    | CONVERSION_STATUS
    | PROCESSING_STATUS;
  type?: FAILED_RETURN_TYPE;
};

export function checkIfFailed(datasetStatus: DatasetUploadStatus): FailReturn {
  if (getConversionStatus(datasetStatus) === CONVERSION_STATUS.FAILED) {
    return {
      error: CONVERSION_STATUS.FAILED,
      isFailed: true,
      type: FAILED_RETURN_TYPE.CONVERSION,
    };
  }

  if (datasetStatus.validation_status === VALIDATION_STATUS.INVALID) {
    return {
      error: VALIDATION_STATUS.INVALID,
      isFailed: true,
      type: FAILED_RETURN_TYPE.VALIDATION,
    };
  }

  if (datasetStatus.upload_status === UPLOAD_STATUS.FAILED) {
    return {
      error: UPLOAD_STATUS.FAILED,
      isFailed: true,
      type: FAILED_RETURN_TYPE.UPLOAD,
    };
  }

  // (thuang): Only show the catch all error state if we don't have better
  // error message to show
  if (datasetStatus.processing_status === PROCESSING_STATUS.FAILURE) {
    return {
      error: PROCESSING_STATUS.FAILURE,
      isFailed: true,
      type: FAILED_RETURN_TYPE.PROCESS,
    };
  }

  return { isFailed: false };
}

export const checkIfLoading = (datasetStatus: DatasetUploadStatus): boolean => {
  if (checkIfFailed(datasetStatus).isFailed) return false;

  if (
    datasetStatus.upload_status === UPLOAD_STATUS.UPLOADING ||
    datasetStatus.upload_status === UPLOAD_STATUS.WAITING
  ) {
    return true;
  }

  if (datasetStatus.validation_status === VALIDATION_STATUS.VALIDATING) {
    return true;
  }

  if (getConversionStatus(datasetStatus) === CONVERSION_STATUS.CONVERTING) {
    return true;
  }

  return datasetStatus.processing_status === PROCESSING_STATUS.PENDING;
};

export function checkIfCancelled(datasetStatus: DatasetUploadStatus): boolean {
  return (
    datasetStatus.upload_status === UPLOAD_STATUS.CANCEL_PENDING ||
    datasetStatus.upload_status === UPLOAD_STATUS.CANCELED
  );
}

export function useCheckCollectionPopulated({
  invalidateCollectionQuery,
  isFailed,
  isNamePopulated,
  validationStatus,
}: {
  invalidateCollectionQuery: () => void;
  isFailed: boolean;
  isNamePopulated: boolean;
  validationStatus: VALIDATION_STATUS;
}): void {
  useCheckCollection({
    invalidateCollectionQuery,
    shouldFetch:
      !isFailed &&
      validationStatus === VALIDATION_STATUS.VALID &&
      !isNamePopulated,
  });
}

type FORMAT_KEYS = "h5ad_status" | "cxg_status" | "rds_status";

const CONVERSION_STATUS_FORMAT_KEYS = [
  "h5ad_status",
  "cxg_status",
  "rds_status",
] as FORMAT_KEYS[];

export function useCheckCollectionFormatsPopulated({
  invalidateCollectionQuery,
  dataset,
  datasetUploadStatus,
}: {
  invalidateCollectionQuery: () => void;
  dataset: Dataset;
  datasetUploadStatus: DatasetUploadStatus;
}): void {
  let numOfConvertedFormats = 0;

  for (const key of CONVERSION_STATUS_FORMAT_KEYS) {
    if (datasetUploadStatus[key] !== CONVERSION_STATUS.CONVERTED) continue;

    numOfConvertedFormats += 1;
  }

  const numOfAvailableFormats = dataset.dataset_assets.length;

  useCheckCollection({
    invalidateCollectionQuery,
    shouldFetch: numOfConvertedFormats > numOfAvailableFormats,
  });
}

const FETCH_COLLECTION_INTERVAL_MS = 5 * 1000;

function useCheckCollection({
  invalidateCollectionQuery,
  shouldFetch,
}: {
  invalidateCollectionQuery: () => void;
  shouldFetch: boolean;
}) {
  useEffect(() => {
    let intervalId: number | undefined = undefined;

    if (!intervalId && shouldFetch) {
      intervalId = window?.setInterval(() => {
        invalidateCollectionQuery();
      }, FETCH_COLLECTION_INTERVAL_MS);
    }

    return () => {
      clearInterval(intervalId);
    };
  }, [invalidateCollectionQuery, shouldFetch]);
}

export function useCancelDatasetStatusQuery({
  datasetId,
  isFailed,
  isLoading,
  queryClient,
}: {
  datasetId: string;
  isFailed: boolean;
  isLoading: boolean;
  queryClient: QueryClient;
}): void {
  useEffect(() => {
    if (isFailed || !isLoading) {
      queryClient.cancelQueries([USE_DATASET_STATUS, datasetId]);
    }
  }, [datasetId, isFailed, isLoading, queryClient]);
}

export function useUploadProgress({
  initProgress,
  progress,
}: {
  initProgress: number;
  progress: number;
}): void {
  useProgress({
    initProgress,
    progress,
    successCallback: useUploadProgressSuccessCallback,
  });
}

function useUploadProgressSuccessCallback() {
  successToast(
    "Upload was successful. Your file is being processed which will " +
      "continue in the background, even if you close this window."
  );
}

export function getConversionStatus(
  datasetStatus?: DatasetUploadStatus
): CONVERSION_STATUS {
  if (!datasetStatus) return CONVERSION_STATUS.NA;

  const { h5ad_status, cxg_status, rds_status } = datasetStatus;

  const statuses = [h5ad_status, cxg_status, rds_status];

  if (statuses.some((status) => status === CONVERSION_STATUS.CONVERTING)) {
    return CONVERSION_STATUS.CONVERTING;
  }

  if (statuses.some((status) => status === CONVERSION_STATUS.FAILED)) {
    return CONVERSION_STATUS.FAILED;
  }

  if (
    statuses.every((status) =>
      [CONVERSION_STATUS.CONVERTED, CONVERSION_STATUS.SKIPPED].includes(status)
    )
  ) {
    return CONVERSION_STATUS.CONVERTED;
  }

  return CONVERSION_STATUS.NA;
}

export function useConversionProgress({
  initDatasetStatus,
  datasetStatus,
}: {
  initDatasetStatus: DatasetUploadStatus;
  datasetStatus?: DatasetUploadStatus;
}): void {
  const initProgress = isConverted(initDatasetStatus);
  const progress = isConverted(datasetStatus);

  useProgress({
    initProgress,
    progress,
    successCallback: useConversionProgressSuccessCallback,
  });
}

function useConversionProgressSuccessCallback() {
  successToast("Processing complete.");
}

function isConverted(status?: DatasetUploadStatus): 0 | 1 {
  return getConversionStatus(status) === CONVERSION_STATUS.CONVERTED ? 1 : 0;
}

export function hasCXGFile(dataset: Dataset): boolean {
  const deployments = dataset.dataset_deployments;

  return !(!deployments || !deployments.length);
}

function successToast(message: string) {
  Toast.show({
    icon: IconNames.TICK,
    intent: Intent.SUCCESS,
    message,
  });
}

function useProgress({
  initProgress,
  progress,
  successCallback,
}: {
  initProgress: number;
  progress: number;
  successCallback: () => void;
}) {
  useEffect(() => {
    if (initProgress === 1) return;

    if (progress === 1) {
      successCallback();
    }
  }, [successCallback, initProgress, progress]);
}
