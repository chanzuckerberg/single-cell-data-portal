import { Button, Intent } from "@blueprintjs/core";
import React, { FC } from "react";
import DeleteDataset from "../DeleteDataset";
import { DatasetStatusTag, StatusContainer, StyledSpinner } from "./style";

interface Props {
  isConverting: boolean;
  isValidating: boolean;
  progress: number;
  datasetId: string;
  collectionId: string;
}

const UploadStatus: FC<Props> = ({
  isConverting,
  isValidating,
  progress,
  datasetId,
  collectionId,
}) => {
  let content: { progress?: number; text: string } = {
    progress,
    text: `Uploading (${Math.round(progress * 100)}%)`,
  };

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
      <DeleteDataset
        id={datasetId}
        collectionId={collectionId}
        Button={CancelButton}
      />
    </StatusContainer>
  );
};

function CancelButton({ ...props }) {
  return (
    <Button minimal {...props}>
      Cancel
    </Button>
  );
}

export default UploadStatus;
