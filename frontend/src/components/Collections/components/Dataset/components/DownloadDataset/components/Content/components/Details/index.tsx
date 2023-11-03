import { Spinner, SpinnerSize } from "@blueprintjs/core";
import { FC, ReactNode } from "react";
import { DialogLoader } from "src/components/Datasets/components/DownloadDataset/style";
import { FormControl, FormLabel } from "@mui/material";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { FEATURES } from "src/common/featureFlags/features";

export const PROMPT_TEXT =
  "Select one of the data formats to view its download details.";

interface Props {
  downloadPreview?: ReactNode;
  selected: boolean;
  fileSize: number;
  isLoading: boolean;
}

const MEGA_BYTES = 2 ** 20;

const Details: FC<Props> = ({
  downloadPreview,
  selected = false,
  fileSize = 0,
  isLoading = false,
}) => {
  const isDownloadUX = useFeatureFlag(FEATURES.DOWNLOAD_UX);
  const Loader = isDownloadUX ? (
    <DialogLoader sdsStyle="minimal" />
  ) : (
    <Spinner size={SpinnerSize.SMALL} />
  ); // TODO: #5566 hidden under feature flag.

  function renderContent() {
    if (isLoading) {
      return Loader;
    }

    if (!selected) {
      return <div>{PROMPT_TEXT}</div>;
    }

    return <div>{`${Math.round(fileSize / MEGA_BYTES)}MB`}</div>;
  }

  return (
    <FormControl>
      <FormLabel>Download Details</FormLabel>
      {renderContent()}
      {downloadPreview}
    </FormControl>
  );
};

export default Details;
