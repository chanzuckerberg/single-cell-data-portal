import React, { FC } from "react";
import Modal from "src/components/common/Modal";
import { SmallColumn } from "../../common/style";
import Content from "./components/Content";
import { StyledButton } from "./style";

const DownloadDataset: FC = () => {
  const [isOpen, setIsOpen] = React.useState(false);

  const toggleOpen = () => setIsOpen(!isOpen);

  return (
    <SmallColumn>
      <StyledButton data-test-id="dataset-download-button" onClick={toggleOpen}>
        Download{" "}
      </StyledButton>
      <Modal isOpen={isOpen} onClose={toggleOpen}>
        <Content />
      </Modal>
    </SmallColumn>
  );
};

export default DownloadDataset;
