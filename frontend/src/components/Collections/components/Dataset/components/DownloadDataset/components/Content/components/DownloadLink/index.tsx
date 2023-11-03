import { FC } from "react";
import { CodeBlock as DownloadUXCodeBlock, DownloadCodeBlock } from "./style";
import CopyButton from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/DownloadLink/components/CopyButton";
import CopyMask from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/DownloadLink/components/CopyMask";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { FEATURES } from "src/common/featureFlags/features";
import { DATASET_ASSET_FORMAT } from "src/common/entities";
import CopyCaption from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/DownloadLink/components/CopyCaption";

interface Props {
  downloadLink: string;
  fileName: string;
  handleAnalytics: () => void;
  selectedFormat: DATASET_ASSET_FORMAT | "";
}

const DownloadLink: FC<Props> = ({
  downloadLink,
  fileName,
  handleAnalytics,
  selectedFormat,
}) => {
  const isDownloadUX = useFeatureFlag(FEATURES.DOWNLOAD_UX);
  const CodeBlock = isDownloadUX ? DownloadUXCodeBlock : DownloadCodeBlock; // TODO(cc) Download UI #5566 hidden under feature flag.
  const curl = `curl -o ${fileName} "${downloadLink}"`;
  return (
    <>
      <CodeBlock>
        {isDownloadUX ? <code>{downloadLink}</code> : <code>{curl}</code>}
        {isDownloadUX ? (
          <CopyButton
            downloadLink={downloadLink}
            handleAnalytics={handleAnalytics}
          />
        ) : (
          <CopyMask curl={curl} handleAnalytics={handleAnalytics} />
        )}
      </CodeBlock>
      <CopyCaption selectedFormat={selectedFormat} />
    </>
  );
};

export default DownloadLink;
