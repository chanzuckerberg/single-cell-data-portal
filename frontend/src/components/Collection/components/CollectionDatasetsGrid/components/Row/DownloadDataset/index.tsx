import React, { ElementType, FC, useState } from "react";
import { useDatasetAssets } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DownloadDataset/util";
import Content from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content";
import {
  Dialog as DownloadUXDialog,
  DownloadDialog,
} from "src/components/Datasets/components/DownloadDataset/style";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { FEATURES } from "src/common/featureFlags/features";

interface Props {
  Button: ElementType;
  datasetId: string;
  isDisabled?: boolean;
  name: string;
}

/**
 * Fetch dataset assets on click of download button. Based on Dataset/components/DownloadDataset/components/index but
 * fetches dataset assets rather than accepts dataset assets as props. Required for core datasets table as dataset
 * assets are not included in the "light" datasets/index endpoint that backs the table and are therefore fetched on
 * demand.
 */
const DownloadDataset: FC<Props> = ({
  Button,
  datasetId,
  isDisabled = false,
  name,
}) => {
  const isDownloadUX = useFeatureFlag(FEATURES.DOWNLOAD_UX);
  const [isOpen, setIsOpen] = useState<boolean>(false);
  const Dialog = isDownloadUX ? DownloadUXDialog : DownloadDialog; // TODO(cc) Download UI #5566 hidden under feature flag.

  // Fetch the dataset assets on open of download modal.
  const { datasetAssets, isError, isLoading } = useDatasetAssets(
    datasetId,
    isOpen
  );

  return (
    <>
      <Button
        datasetName={name}
        data-testid="dataset-download-button"
        disabled={isDisabled}
        onClick={() => setIsOpen(true)}
      />
      <Dialog onClose={() => setIsOpen(false)} open={isOpen}>
        <Content
          dataAssets={datasetAssets}
          isError={isError}
          isLoading={isLoading}
          name={name}
          onClose={() => setIsOpen(false)}
        />
      </Dialog>
    </>
  );
};

export default DownloadDataset;
