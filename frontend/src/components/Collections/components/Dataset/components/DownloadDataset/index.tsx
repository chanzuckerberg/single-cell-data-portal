import React, { ElementType, FC } from "react";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { Dataset } from "src/common/entities";
import Content from "./components/Content";
import {
  Dialog as DownloadUXDialog,
  DownloadDialog,
} from "src/components/Datasets/components/DownloadDataset/style";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { FEATURES } from "src/common/featureFlags/features";

interface Props {
  Button: ElementType;
  dataAssets: Dataset["dataset_assets"];
  isDisabled?: boolean;
  name: string;
}

const DownloadDataset: FC<Props> = ({
  Button,
  dataAssets = EMPTY_ARRAY,
  isDisabled = false,
  name,
}) => {
  const isDownloadUX = useFeatureFlag(FEATURES.DOWNLOAD_UX);
  const [isOpen, setIsOpen] = React.useState(false);
  const Dialog = isDownloadUX ? DownloadUXDialog : DownloadDialog; // TODO(cc) Download UI #5566 hidden under feature flag.

  if (!dataAssets.length) {
    return null;
  }

  return (
    <>
      <Button
        datasetName={name}
        data-testid="dataset-download-button"
        disabled={isDisabled || !dataAssets.length}
        onClick={() => setIsOpen(true)}
      />
      <Dialog onClose={() => setIsOpen(false)} open={isOpen}>
        <Content
          name={name}
          dataAssets={dataAssets}
          onClose={() => setIsOpen(false)}
        />
      </Dialog>
    </>
  );
};

export default DownloadDataset;
