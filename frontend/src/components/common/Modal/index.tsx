import { Classes } from "@blueprintjs/core";
import React, { FC } from "react";
import { StyledDialog } from "./style";

interface Props {
  onClose: () => void;
  isOpen: boolean;
  title: string;
}

const Modal: FC<Props> = ({ onClose, title, isOpen, children }) => {
  return (
    <StyledDialog title={title} isOpen={isOpen} onClose={onClose}>
      <div className={Classes.DIALOG_BODY}>{children}</div>
    </StyledDialog>
  );
};

export default Modal;
