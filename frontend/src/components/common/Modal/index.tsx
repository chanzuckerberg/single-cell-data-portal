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
      {children}
    </StyledDialog>
  );
};

export default Modal;
