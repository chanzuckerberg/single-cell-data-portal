import React, { FC } from "react";
import Modal from "src/components/common/Modal";
import Content from "./components/Content";

const DownloadDataset: FC = () => {
  const [isOpen, setIsOpen] = React.useState(false);

  const toggleOpen = () => setIsOpen(!isOpen);

  return (
    <div>
      <button onClick={toggleOpen}>Download</button>
      <Modal isOpen={isOpen} onClose={toggleOpen}>
        <Content />
      </Modal>
    </div>
  );
};

export default DownloadDataset;
