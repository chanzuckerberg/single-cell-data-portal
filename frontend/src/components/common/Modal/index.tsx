import { FC } from "react";
import { StyledDialog } from "./style";

interface Props {
  className?: string;
  onClose: () => void;
  isOpen: boolean;
  title: string;
}

const Modal: FC<Props> = ({ className, onClose, title, isOpen, children }) => {
  return (
    <StyledDialog
      title={title}
      isOpen={isOpen}
      onClose={onClose}
      className={className}
      shouldReturnFocusOnClose={false}
    >
      {children}
    </StyledDialog>
  );
};

export default Modal;
