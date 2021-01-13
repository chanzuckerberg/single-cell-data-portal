import { Intent, Spinner } from "@blueprintjs/core";
import React, { FC } from "react";
import { DatasetStatusTag } from "./style";

interface Props {
  isValidating: boolean;
  progress: number;
}

const UploadStatus: FC<Props> = ({ isValidating, progress }) => {
  if (isValidating) {
    return (
      <DatasetStatusTag>
        <Spinner intent={Intent.PRIMARY} size={16} />
        Validating...
      </DatasetStatusTag>
    );
  }

  return (
    <DatasetStatusTag>
      <Spinner intent={Intent.PRIMARY} value={progress} size={16} />
      {`Uploading (${Math.round(progress * 100)}%)`}
    </DatasetStatusTag>
  );
};

export default UploadStatus;
