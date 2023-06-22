import * as React from "react";
import { FC } from "react";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { Dataset } from "src/common/entities";
import Modal from "src/components/common/Modal";
import Content from "./components/Content";
import { StyledButton } from "./style";

interface Props {
  name: string;
  dataAssets: Dataset["dataset_assets"];
  isDisabled?: boolean;
  Button?: React.ElementType;
}

const DownloadDataset: FC<Props> = ({
  name,
  dataAssets = EMPTY_ARRAY,
  isDisabled = false,
  Button = StyledButton,
}) => {
  const [isOpen, setIsOpen] = React.useState(false);

  const toggleOpen = () => {
    setIsOpen(!isOpen);
  };

  if (!dataAssets.length) {
    return null;
  }

  return (
    <>
      <Button
        datasetName={name}
        disabled={isDisabled || !dataAssets.length}
        onClick={toggleOpen}
        data-testid="dataset-download-button"
      />
      <Modal
        isCloseButtonShown={false}
        isOpen={isOpen}
        onClose={toggleOpen}
        title="Download Dataset"
      >
        <Content name={name} dataAssets={dataAssets} onClose={toggleOpen} />
      </Modal>
    </>
  );
};

export default DownloadDataset;
