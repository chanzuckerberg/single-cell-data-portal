import { Intent, Spinner } from "@blueprintjs/core";
import React, { FC } from "react";
import { CancelButton, DatasetStatusTag, StatusContainer } from "./style";

interface Props {
  isValidating: boolean;
  progress: number;
}

const UploadStatus: FC<Props> = ({ isValidating, progress }) => {
  if (isValidating) {
    return (
      <StatusContainer>
        <DatasetStatusTag>
          <Spinner intent={Intent.PRIMARY} size={16} />
          Validating...
        </DatasetStatusTag>
        <CancelButton>Cancel</CancelButton>
      </StatusContainer>
    );
  }

  return (
    <StatusContainer>
      <DatasetStatusTag>
        <Spinner intent={Intent.PRIMARY} value={progress} size={16} />
        {`Uploading (${Math.round(progress * 100)}%)`}
      </DatasetStatusTag>
      <CancelButton>Cancel</CancelButton>
    </StatusContainer>
  );
};

export default UploadStatus;
