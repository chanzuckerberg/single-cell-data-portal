import { Intent, Spinner } from "@blueprintjs/core";
import React, { FC } from "react";
import { CancelButton, DatasetStatusTag, StatusContainer } from "./style";

interface Props {
  isValidating: boolean;
  progress: number;
  cancelUpload: () => void;
}

const UploadStatus: FC<Props> = ({ isValidating, progress, cancelUpload }) => {
  const cancel = <CancelButton onClick={cancelUpload}>Cancel</CancelButton>;

  if (isValidating) {
    return (
      <StatusContainer>
        <DatasetStatusTag>
          <Spinner intent={Intent.PRIMARY} size={16} />
          Validating...
        </DatasetStatusTag>
        {cancel}
      </StatusContainer>
    );
  }

  return (
    <StatusContainer>
      <DatasetStatusTag>
        <Spinner intent={Intent.PRIMARY} value={progress} size={16} />
        {`Uploading (${Math.round(progress * 100)}%)`}
      </DatasetStatusTag>
      {cancel}
    </StatusContainer>
  );
};

export default UploadStatus;
