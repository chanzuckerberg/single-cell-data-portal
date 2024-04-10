import { FC } from "react";
import { CodeBlock } from "./style";
import CopyButton from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/DownloadLink/components/CopyButton";
import { DATASET_ASSET_FORMAT } from "src/common/entities";
import CopyCaption from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/DownloadLink/components/CopyCaption";

interface Props {
  downloadLink: string;
  handleAnalytics: () => void;
  selectedFormat: DATASET_ASSET_FORMAT | "";
}

const DownloadLink: FC<Props> = ({
  downloadLink,
  handleAnalytics,
  selectedFormat,
}) => {
  return (
    <>
      <CodeBlock>
        <code>{downloadLink}</code>
        <CopyButton
          downloadLink={downloadLink}
          handleAnalytics={handleAnalytics}
        />
      </CodeBlock>
      <CopyCaption selectedFormat={selectedFormat} />
    </>
  );
};

export default DownloadLink;
