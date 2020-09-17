import copy from "clipboard-copy";
import React, { FC, useState } from "react";
import { Code, CodeMask, CodeWrapper, Tip } from "./style";

interface Props {
  link: string;
  fileName: string;
}

const CurlLink: FC<Props> = ({ link, fileName }) => {
  const [isCopied, setIsCopied] = useState(false);

  const curl = `curl -o ${fileName} "${link}"`;

  const handleCopyClick = () => {
    setIsCopied(true);
    copy(curl);
  };
  const handleCopyMouseEnter = () => setIsCopied(false);

  return (
    <>
      <CodeWrapper>
        <Code>{curl}</Code>
        <CodeMask onClick={handleCopyClick} onMouseEnter={handleCopyMouseEnter}>
          {isCopied ? "Copied!" : "Copy to Clipboard"}
        </CodeMask>
      </CodeWrapper>
      <Tip>
        If you prefer not to download this dataset directly in your browser, you
        can optionally use the provided cURL link to download via the terminal.
        The above link will be valid for 1 week.
      </Tip>
    </>
  );
};

export default CurlLink;
