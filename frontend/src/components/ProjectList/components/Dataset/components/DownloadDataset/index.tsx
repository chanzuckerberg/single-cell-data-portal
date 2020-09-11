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
        <StyledButton onClick={toggleOpen}>
          Download
          <Modal isOpen={isOpen} onClose={toggleOpen}>
            <Content />
          </Modal>
        </StyledButton>
      </Wrapper>
    </div>
  );
};

export default DownloadDataset;
