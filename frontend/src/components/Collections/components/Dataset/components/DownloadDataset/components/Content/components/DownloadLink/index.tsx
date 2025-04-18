import { FC } from "react";
import { CodeBlock } from "./style";
import CopyButton from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/DownloadLink/components/CopyButton";
import CopyCaption from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/DownloadLink/components/CopyCaption";
import type { DownloadLinkType } from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content";
interface Props {
  downloadLinks: DownloadLinkType[];
  handleAnalytics: () => void;
  formatsToDownload: string[];
}

const DownloadLink: FC<Props> = ({
  downloadLinks,
  handleAnalytics,
  formatsToDownload,
}) => {
  const copyText = downloadLinks.reduce((acc, download) => {
    if (!formatsToDownload.includes(download.filetype)) {
      return acc;
    }
    return acc + download.downloadURL + "\n";
  }, "");

  return (
    <>
      <CodeBlock>
        <code>{copyText}</code>
        <CopyButton
          downloadLink={copyText}
          handleAnalytics={handleAnalytics}
          label={downloadLinks.length > 1 ? "Copy All" : "Copy"}
        />
      </CodeBlock>
      <CopyCaption />
    </>
  );
};

export default DownloadLink;
