import React, { ElementType, FC } from "react";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { Dataset } from "src/common/entities";
import Content from "./components/Content";
import { Dialog } from "src/components/Datasets/components/DownloadDataset/style";
import { useDialog } from "src/views/Collection/hooks/useDialog";

interface Props {
  Button: ElementType;
  dataAssets: Dataset["dataset_assets"];
  isDisabled?: boolean;
  name: string;
}

const DownloadDataset: FC<Props> = ({
  Button,
  dataAssets = EMPTY_ARRAY,
  isDisabled = false,
  name,
}) => {
  const { onClose, onOpen, open } = useDialog();

  if (!dataAssets.length) {
    return null;
  }

  return (
    <>
      <Button
        datasetName={name}
        data-testid="dataset-download-button"
        disabled={isDisabled || !dataAssets.length}
        onClick={onOpen}
      />
      <Dialog onClose={onClose} open={open}>
        <Content name={name} dataAssets={dataAssets} onClose={onClose} />
      </Dialog>
    </>
  );
};

export default DownloadDataset;
