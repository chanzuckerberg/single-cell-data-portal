import React, { FC } from "react";
import { Dataset } from "src/common/entities";
import Modal from "src/components/common/Modal";
import { SmallColumn } from "../../common/style";
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
  dataAssets,
  isDisabled = false,
  Button = StyledButton,
}) => {
  const [isOpen, setIsOpen] = React.useState(false);

  const toggleOpen = () => setIsOpen(!isOpen);

  return (
    <SmallColumn>
      <Button
        disabled={isDisabled}
        onClick={toggleOpen}
        data-test-id="dataset-download-button"
      >
        Download
      </Button>
      <Modal title="Download Dataset" isOpen={isOpen} onClose={toggleOpen}>
        <Content name={name} dataAssets={dataAssets} onClose={toggleOpen} />
      </Modal>
    </SmallColumn>
  );
};

export default DownloadDataset;
