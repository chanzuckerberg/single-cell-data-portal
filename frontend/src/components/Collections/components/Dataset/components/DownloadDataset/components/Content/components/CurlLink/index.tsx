import { Link } from "@czi-sds/components";
import copy from "clipboard-copy";
import { FC, useState } from "react";
import { Code, CodeMask, CodeWrapper, Tip } from "./style";

const DISCOVER_API_URL = "https://api.cellxgene.cziscience.com/curation/ui/#/";
const SCHEMA_URL =
  "https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.1.0/schema.md";

interface Props {
  fileName: string;
  handleAnalytics: () => void;
  link: string;
}

const CurlLink: FC<Props> = ({ fileName, handleAnalytics, link }) => {
  const [isCopied, setIsCopied] = useState(false);

  const curl = `curl -o ${fileName} "${link}"`;

  const handleCopyClick = () => {
    setIsCopied(true);
    copy(curl);
    handleAnalytics();
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
        This cURL command is valid forever. All datasets on CZ CELLxGENE
        Discover adhere to its{" "}
        <Link href={SCHEMA_URL} rel="noreferrer noopener" target="_blank">
          schema
        </Link>{" "}
        and may be downloaded programmatically using the{" "}
        <Link href={DISCOVER_API_URL} rel="noreferrer noopener" target="_blank">
          Discover API
        </Link>
        .
      </Tip>
    </>
  );
};

export default CurlLink;
