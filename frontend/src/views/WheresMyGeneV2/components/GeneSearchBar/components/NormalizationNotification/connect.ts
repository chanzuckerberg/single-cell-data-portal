import { useEffect, useState } from "react";

export const useConnect = () => {
  const [showURLCopyNotification, setShowURLCopyNotification] = useState(0);

  useEffect(() => {
    if (showURLCopyNotification) {
      setShowURLCopyNotification((prev) => prev + 1);
    }
  }, [showURLCopyNotification]);
  return { showURLCopyNotification };
};
