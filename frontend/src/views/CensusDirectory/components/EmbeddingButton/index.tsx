import React from "react";
import { StyledButton, StyledRadioGroup } from "../../styles";

import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  InputRadio,
  DialogTitle,
} from "@czi-sds/components";
import { useConnect } from "./connect";
import { CopyMask } from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/DownloadLink/components/CopyMask/style";
import { DownloadCodeBlock } from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/DownloadLink/style";

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
          <DownloadCodeBlock>
            <code>{language === "python" ? pythonCode : rCode}</code>
            <CopyMask
              onClick={handleCopyClick}
              onMouseEnter={handleCopyMouseEnter}
            >
              {isCopied ? "Copied!" : "Copy to Clipboard"}
            </CopyMask>
          </DownloadCodeBlock>
        </DialogContent>
        <DialogActions>
          <Button onClick={handleButtonClick}>Close</Button>
        </DialogActions>
      </Dialog>
    </>
  );
}

export default EmbeddingButton;
