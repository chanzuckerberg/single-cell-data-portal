import { Dialog } from "@blueprintjs/core";
import React, { FC } from "react";

interface Props {
  onClose: () => void;
  isOpen: boolean;
}

const Modal: FC<Props> = ({ onClose, isOpen, children }) => {
  return (
    <div>
      <Dialog isOpen={isOpen} onClose={onClose}>
        <div>{children}</div>
      </Dialog>
    </div>
  );
};

export default Modal;
