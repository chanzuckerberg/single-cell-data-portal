import { Intent } from "@blueprintjs/core";
import React, { FC } from "react";
import { DatasetStatusTag, StatusContainer, StyledSpinner } from "./style";

interface Props {
  isConverting: boolean;
  isValidating: boolean;
  isWaiting: boolean;
  progress: number;
}

const UploadStatus: FC<Props> = ({
  isConverting,
  isValidating,
  isWaiting,
  progress,
}) => {
  let content: { progress?: number; text: string } = {
    progress,
    text: `Uploading (${Math.round(progress * 100)}%)`,
  };

  if (isWaiting) {
    content = {
      text: "Queued...",
    };
  }

  if (isConverting) {
    content = {
      text: "Processing...",
    };
  }

  if (isValidating) {
    content = {
      text: "Validating...",
    };
  }

  return (
    <StatusContainer>
      <DatasetStatusTag>
        <StyledSpinner
          intent={Intent.PRIMARY}
          value={content.progress}
          size={16}
        />
        {content.text}
      </DatasetStatusTag>
    </StatusContainer>
  );
};

export default UploadStatus;
