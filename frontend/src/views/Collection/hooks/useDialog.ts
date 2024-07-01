import { useCallback, useState } from "react";

export interface UseDialog {
  onClose: () => void;
  onOpen: () => void;
  open: boolean;
}

/**
 * Dialog state (open / close) functionality for collection view dialog components.
 * @returns dialog state functionality.
 */
export function useDialog(): UseDialog {
  const [open, setOpen] = useState<boolean>(false);

  // Opens dialog.
  const onOpen = useCallback(() => {
    setOpen(true);
  }, []);

  // Closes dialog.
  const onClose = useCallback(() => {
    setOpen(false);
  }, []);

  return { onClose, onOpen, open };
}
