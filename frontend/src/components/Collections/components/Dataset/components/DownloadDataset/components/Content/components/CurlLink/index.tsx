import { FC } from "react";
import {
  Caption as DownloadUXCaption,
  CodeBlock as DownloadUXCodeBlock,
  DownloadCaption,
  DownloadCodeBlock,
} from "./style";
import CopyButton from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/CurlLink/components/CopyButton";
import CopyMask from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/CurlLink/components/CopyMask";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { FEATURES } from "src/common/featureFlags/features";

interface Props {
  fileName: string;
  handleAnalytics: () => void;
  link: string;
}

const CurlLink: FC<Props> = ({ fileName, handleAnalytics, link }) => {
  const isDownloadUX = useFeatureFlag(FEATURES.DOWNLOAD_UX);
  const CodeBlock = isDownloadUX ? DownloadUXCodeBlock : DownloadCodeBlock; // TODO(cc) Download UI #5566 hidden under feature flag.
  const Caption = isDownloadUX ? DownloadUXCaption : DownloadCaption; // TODO(cc) Download UI #5566 hidden under feature flag.
  const curl = `curl -o ${fileName} "${link}"`;
  return (
    <>
      <CodeBlock>
        <code>{curl}</code>
        {isDownloadUX ? (
          <CopyButton curl={curl} handleAnalytics={handleAnalytics} />
        ) : (
          <CopyMask curl={curl} handleAnalytics={handleAnalytics} />
        )}
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
