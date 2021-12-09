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
      // Remove this prop after upgrading Blueprint Tooltip library
      // Current bug is after closing the modal, tooltip will reappear, as shown in the link below:
      // https://codesandbox.io/s/blueprint-tooltip-bug-forked-jw9z2?file=/src/index.tsx
      shouldReturnFocusOnClose={false}
    >
      {children}
    </StyledDialog>
  );
};

export default Modal;
