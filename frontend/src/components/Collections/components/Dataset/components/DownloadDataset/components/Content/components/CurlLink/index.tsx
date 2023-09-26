import { FC } from "react";
import { Caption, CodeBlock } from "./style";
import CopyButton from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/CurlLink/components/CopyButton";

interface Props {
  fileName: string;
  handleAnalytics: () => void;
  link: string;
}

const CurlLink: FC<Props> = ({ fileName, handleAnalytics, link }) => {
  const curl = `curl -o ${fileName} "${link}"`;
  return (
    <>
      <CodeBlock>
        <code>{curl}</code>
        <CopyButton curl={curl} handleAnalytics={handleAnalytics} />
      </CodeBlock>
      <Caption>
        If you prefer not to download this dataset directly in your browser, you
        can optionally use the provided cURL link to download via the terminal.
        The above link will be valid for 1 week.
      </Caption>
    </>
  );
};

export default CurlLink;
