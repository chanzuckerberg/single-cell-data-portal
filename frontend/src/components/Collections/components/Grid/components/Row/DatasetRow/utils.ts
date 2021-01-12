import {
  CONVERSION_STATUS,
  DatasetUploadStatus,
  UPLOAD_STATUS,
  VALIDATION_STATUS,
} from "src/common/entities";

export type FailReturn = {
  isFailed: boolean;
  error?: VALIDATION_STATUS | UPLOAD_STATUS;
};

export function checkIfFailed(datasetStatus: DatasetUploadStatus): FailReturn {
  if (datasetStatus.validation_status === VALIDATION_STATUS.INVALID) {
    return { error: VALIDATION_STATUS.INVALID, isFailed: true };
  }

  if (datasetStatus.upload_status === UPLOAD_STATUS.FAILED) {
    return { error: UPLOAD_STATUS.FAILED, isFailed: true };
  }

  // TODO: check if conversion failed
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

  if (
    datasetStatus.conversion_anndata_status === CONVERSION_STATUS.CONVERTING ||
    datasetStatus.conversion_cxg_status === CONVERSION_STATUS.CONVERTING ||
    datasetStatus.conversion_rds_status === CONVERSION_STATUS.CONVERTING ||
    datasetStatus.conversion_loom_status === CONVERSION_STATUS.CONVERTING
  ) {
    return true;
  }

  return false;
};
