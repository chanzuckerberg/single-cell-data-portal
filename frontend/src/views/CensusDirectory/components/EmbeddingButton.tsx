import React, { useCallback, useState } from "react";
import { StyledButton, StyledRadioGroup } from "../styles";
import {
  Code,
  CodeMask,
  CodeWrapper,
} from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/CurlLink/style";
import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  InputRadio,
  DialogTitle,
} from "@czi-sds/components";

function EmbeddingButton() {
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

  return (
    <>
      <StyledButton
        sdsType="secondary"
        sdsStyle="square"
        onClick={handleButtonClick}
      >
        Embedding
      </StyledButton>
      <Dialog open={isOpen} onClose={handleButtonClick}>
        <DialogTitle title="Embedding" />
        <DialogContent>
          To get started with this embedding, paste the following into your
          code:
          <StyledRadioGroup
            value={language}
            onChange={(_, value: string) => setLanguage(value)}
            row
          >
            <InputRadio label="Python" value="python" />
            <InputRadio label="R" value="r" />
          </StyledRadioGroup>
          <CodeWrapper>
            <Code>{language === "python" ? pythonCode : rCode}</Code>
            <CodeMask
              onClick={handleCopyClick}
              onMouseEnter={handleCopyMouseEnter}
            >
              {isCopied ? "Copied!" : "Copy to Clipboard"}
            </CodeMask>
          </CodeWrapper>
        </DialogContent>
        <DialogActions>
          <Button onClick={handleButtonClick}>Close</Button>
        </DialogActions>
      </Dialog>
    </>
  );
}

export default EmbeddingButton;
