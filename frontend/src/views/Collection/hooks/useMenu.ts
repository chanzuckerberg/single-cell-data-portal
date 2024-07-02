import { MouseEvent, useCallback, useState } from "react";

export interface UseMenu {
  anchorEl: HTMLButtonElement | null;
  onClose: () => void;
  onOpen: (mouseEvent: MouseEvent<HTMLButtonElement>) => void;
  open: boolean;
}

/**
 * Menu functionality for collection view menu dropdown.
 * @returns menu functionality.
 */
export function useMenu(): UseMenu {
  const [anchorEl, setAnchorEl] = useState<null | HTMLButtonElement>(null);
  const open = Boolean(anchorEl);

  // Opens menu.
  const onOpen = useCallback((mouseEvent: MouseEvent<HTMLButtonElement>) => {
    setAnchorEl(mouseEvent.currentTarget);
  }, []);

  // Closes menu.
  const onClose = useCallback(() => {
    setAnchorEl(null);
  }, []);

  return { anchorEl, onClose, onOpen, open };
}
