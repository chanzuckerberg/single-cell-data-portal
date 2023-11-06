import React from "react";
import { StyledButton, StyledRadioGroup } from "../../styles";
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
import { useConnect } from "./connect";

function EmbeddingButton() {
  const {
    isOpen,
    isCopied,
    language,
    pythonCode,
    rCode,
    setLanguage,
    handleButtonClick,
    handleCopyClick,
    handleCopyMouseEnter,
  } = useConnect();
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
