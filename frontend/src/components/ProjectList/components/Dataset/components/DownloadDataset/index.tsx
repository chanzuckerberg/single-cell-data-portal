import React, { FC } from "react";
import Modal from "src/components/common/Modal";
import Content from "./components/Content";

const DownloadDataset: FC = () => {
  // DEBUG
  // DEBUG
  // DEBUG
  const [isOpen, setIsOpen] = React.useState(true);

  const toggleOpen = () => setIsOpen(!isOpen);

  return (
    <div>
      <button onClick={toggleOpen}>Download</button>
      <Modal title="Download Dataset" isOpen={isOpen} onClose={toggleOpen}>
        <Content />
      </Modal>
    </div>
  );
};

export default DownloadDataset;
