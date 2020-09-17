import React, { FC } from "react";
import { Dataset } from "src/common/entities";
import Modal from "src/components/common/Modal";
import { SmallColumn } from "../../common/style";
import Content from "./components/Content";
import { StyledButton } from "./style";

interface Props {
  name: string;
  dataAssets: Dataset["dataset_assets"];
}

const DownloadDataset: FC<Props> = ({ name, dataAssets }) => {
  const [isOpen, setIsOpen] = React.useState(false);

  const toggleOpen = () => setIsOpen(!isOpen);

  return (
    <SmallColumn>
      <StyledButton onClick={toggleOpen} data-test-id="dataset-download-button">
        Download
      </StyledButton>
      <Modal title="Download Dataset" isOpen={isOpen} onClose={toggleOpen}>
        <Content name={name} dataAssets={dataAssets} onClose={toggleOpen} />
      </Modal>
    </SmallColumn>
  );
};

export default DownloadDataset;
