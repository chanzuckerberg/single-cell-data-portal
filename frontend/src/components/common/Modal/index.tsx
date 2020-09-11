import { Classes, Overlay } from "@blueprintjs/core";
import React, { FC } from "react";

interface Props {
  onClose: () => void;
  isOpen: boolean;
}

const Modal: FC<Props> = ({ onClose, isOpen, children }) => {
  return (
    <div>
      <Overlay isOpen={isOpen} onClose={onClose}>
        <div className={Classes.OVERLAY_CONTENT}>{children}</div>
      </Overlay>
    </div>
  );
};

export default Modal;
