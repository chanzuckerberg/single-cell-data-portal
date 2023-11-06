import { useCallback, useState } from "react";

export const useConnect = () => {
  const [isOpen, setIsOpen] = useState(false);
  const [isCopied, setIsCopied] = useState(false);
  const [language, setLanguage] = useState<string>("python");

  const handleButtonClick = useCallback(() => {
    // TODO: Analytics
    // if (!isOpen) track(EVENTS.WMG_DOWNLOAD_CLICKED);
    setIsOpen(!isOpen);
  }, [isOpen]);

  // These can be derived from the static S3 namespace + the accessor_id or will be a static url provided in json blob
  const pythonCode = "Long-arbitrary-string-here-python";
  const rCode = "Long-arbitrary-string-here-r";

  const handleCopyClick = () => {
    setIsCopied(true);
    navigator.clipboard.writeText(language === "python" ? pythonCode : rCode);
    // TODO: analytics
  };
  const handleCopyMouseEnter = () => setIsCopied(false);

  return {
    isOpen,
    isCopied,
    language,
    pythonCode,
    rCode,
    setLanguage,
    handleButtonClick,
    handleCopyClick,
    handleCopyMouseEnter,
  };
};
