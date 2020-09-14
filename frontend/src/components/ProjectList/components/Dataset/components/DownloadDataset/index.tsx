import React, { FC } from "react";
import Modal from "src/components/common/Modal";
import Content from "./components/Content";
import { StyledButton, Wrapper } from "./style";

const DownloadDataset: FC = () => {
  const [isOpen, setIsOpen] = React.useState(false);

  const toggleOpen = () => setIsOpen(!isOpen);

  return (
    <div>
      <Wrapper>
        <StyledButton
          data-test-id="dataset-download-button"
          onClick={toggleOpen}
        >
          Download{" "}
        </StyledButton>
        <Modal isOpen={isOpen} onClose={toggleOpen}>
          <Content />
        </Modal>
      </Wrapper>
    </div>
  );
};

export default DownloadDataset;
