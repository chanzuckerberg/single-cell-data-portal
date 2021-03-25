import { Button, Intent } from "@blueprintjs/core";
import React, { FC } from "react";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { Collection, Dataset, VISIBILITY_TYPE } from "src/common/entities";
import { hasAssets } from "src/common/modules/datasets/selectors";
import { useCollection } from "src/common/queries/collections";
import DownloadDataset from "src/components/Collections/components/Dataset/components/DownloadDataset";
import DeleteDataset from "src/components/Collections/components/Grid/components/Row/DatasetRow/components/DeleteDataset";
import DropboxChooser, {
  Props as DropboxChooserProps,
  UploadingFile,
} from "src/components/DropboxChooser";
import { StyledDiv } from "./style";
import { DownloadButton, getSelectedDataset } from "./utils";

export interface UploadedFiles {
  [datasetID: string]: UploadingFile;
}

interface Props {
  collectionId: Collection["id"];
  selectedDatasetId: Dataset["id"];
  visibility: VISIBILITY_TYPE;
  addNewFile: DropboxChooserProps["onUploadFile"];
}

const ActionButtons: FC<Props> = ({
  collectionId,
  selectedDatasetId,
  visibility,
  addNewFile,
}) => {
  const { data: collection } = useCollection({ id: collectionId, visibility });

  const selectedDataset = getSelectedDataset({
    datasets: collection?.datasets,
    selectedId: selectedDatasetId,
  });

  return (
    <StyledDiv>
      {visibility === VISIBILITY_TYPE.PRIVATE && (
        <DropboxChooser onUploadFile={addNewFile}>
          <Button intent={Intent.PRIMARY} outlined>
            Add
          </Button>
        </DropboxChooser>
      )}

      <DownloadDataset
        isDisabled={!hasAssets(selectedDataset)}
        name={selectedDataset?.name || ""}
        dataAssets={selectedDataset?.dataset_assets || EMPTY_ARRAY}
        Button={DownloadButton}
      />

      {visibility === VISIBILITY_TYPE.PRIVATE && (
        <DeleteDataset id={selectedDatasetId} collectionId={collection?.id} />
      )}
    </StyledDiv>
  );
};

export default ActionButtons;
