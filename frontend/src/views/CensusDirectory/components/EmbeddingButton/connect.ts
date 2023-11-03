import { useCallback, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { EmbeddingButtonProps } from "./types";

export const useConnect = ({ project }: EmbeddingButtonProps) => {
  const [isOpen, setIsOpen] = useState(false);
  const [isCopied, setIsCopied] = useState(false);
  const [language, setLanguage] = useState<string>("python");

  const handleButtonClick = useCallback(() => {
    if (!isOpen)
      track(EVENTS.CENSUS_EMBEDDING_CLICKED, {
        project: project.title,
        category: "tier" in project ? project.tier : "2",
      });
    setIsOpen(!isOpen);
  }, [isOpen, project]);

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
