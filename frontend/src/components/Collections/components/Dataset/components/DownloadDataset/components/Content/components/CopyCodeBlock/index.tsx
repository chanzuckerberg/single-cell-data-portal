import { FC, ReactNode } from "react";
import { CodeBlock } from "./style";
import CopyButton from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/CopyCodeBlock/components/CopyButton";
import CopyCaption from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/CopyCodeBlock/components/CopyCaption";
interface Props {
  plural?: boolean;
  handleAnalytics: () => void;
  copyText: string;
  captionText: ReactNode;
}

const CopyCodeBlock: FC<Props> = ({
  plural = false,
  handleAnalytics,
  copyText,
  captionText,
}) => {
  return (
    <>
      <CodeBlock>
        <code>{copyText}</code>
        <CopyButton
          downloadLink={copyText}
          handleAnalytics={handleAnalytics}
          label={plural ? "Copy All" : "Copy"}
        />
      </CodeBlock>
      <CopyCaption>{captionText}</CopyCaption>
    </>
  );
};

export default CopyCodeBlock;
