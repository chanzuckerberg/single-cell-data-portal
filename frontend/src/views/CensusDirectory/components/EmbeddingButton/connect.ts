import { useCallback, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { EmbeddingButtonProps } from "./types";
import { getProjectTier } from "../../utils";

export const useConnect = ({ project }: EmbeddingButtonProps) => {
  const [isOpen, setIsOpen] = useState(false);
  const [isCopied, setIsCopied] = useState(false);
  const [language, setLanguage] = useState<string>("python");

  const projectTier = getProjectTier(project);

  const handleButtonClick = useCallback(() => {
    if (!isOpen)
      track(EVENTS.CENSUS_EMBEDDING_CLICKED, {
        project: project.title,
        category: projectTier,
      });
    setIsOpen(!isOpen);
  }, [isOpen, projectTier, project.title]);

  // These can be derived from the static S3 namespace + the accessor_id or will be a static url provided in json blob
  const pythonCode = "Long-arbitrary-string-here-python";
  const rCode = "Long-arbitrary-string-here-r";

  const handleCopyClick = () => {
    setIsCopied(true);
    navigator.clipboard.writeText(language === "python" ? pythonCode : rCode);
    track(EVENTS.CENSUS_EMBEDDING_COPIED, {
      project: project.title,
      category: projectTier,
      version: language,
    });
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
