import { ReactNode } from "react";
import { StyledDialog } from "./style";

interface Props {
  className?: string;
  onClose: () => void;
  isOpen: boolean;
  title?: string;
  isCloseButtonShown?: boolean;
  children: ReactNode;
}

const Modal = ({
  className,
  onClose,
  title,
  isOpen,
  children,
  isCloseButtonShown,
}: Props) => {
  return (
    <StyledDialog
      title={title}
      isOpen={isOpen}
      onClose={onClose}
      className={className}
      isCloseButtonShown={isCloseButtonShown}
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
