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
import { Link } from "@czi-sds/components";

const DISCOVER_API_URL = "https://api.cellxgene.cziscience.com/curation/ui/#/";
const SCHEMA_URL =
  "https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.1.0/schema.md";

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
      </Caption>
    </>
  );
};

export default CurlLink;
