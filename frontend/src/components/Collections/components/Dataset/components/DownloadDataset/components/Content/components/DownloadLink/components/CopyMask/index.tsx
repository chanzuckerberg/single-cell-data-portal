import copy from "clipboard-copy";
import { useState } from "react";
import { CopyMask as Mask } from "./style";

interface Props {
  curl: string;
  handleAnalytics: () => void;
}

/**
 * @deprecated by CopyButton component once "DOWNLOAD_UX" feature flag is removed (#5566).
 */
export default function CopyMask({
  curl,
  handleAnalytics,
}: Props): JSX.Element {
  const [isCopied, setIsCopied] = useState<boolean>(false);

  // Copy to clipboard, handle analytics.
  const handleCopyClick = () => {
    setIsCopied(true);
    copy(curl);
    handleAnalytics();
  };

  // Handle mouse enter event to reset the copy state.
  const onMouseEnter = () => {
    setIsCopied(false);
  };

  return (
    <Mask onClick={handleCopyClick} onMouseEnter={onMouseEnter}>
      {isCopied ? "Copied!" : "Copy to Clipboard"}
    </Mask>
  );
}
