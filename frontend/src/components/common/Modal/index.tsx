import { Classes } from "@blueprintjs/core";
import React, { FC } from "react";
import { StyledDialog } from "./style";

interface Props {
  onClose: () => void;
  isOpen: boolean;
}

const Modal: FC<Props> = ({ onClose, isOpen, children }) => {
  return (
    <StyledDialog isOpen={isOpen} onClose={onClose}>
      <div className={Classes.DIALOG_BODY}>{children}</div>
    </StyledDialog>
  );
};

export default Modal;
