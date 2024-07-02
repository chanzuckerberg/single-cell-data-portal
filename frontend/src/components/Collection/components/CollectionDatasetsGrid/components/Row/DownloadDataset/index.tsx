import React, { ElementType, FC } from "react";
import { useDatasetAssets } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DownloadDataset/util";
import Content from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content";
import { StyledDialog } from "src/components/Datasets/components/DownloadDataset/style";
import { useDialog } from "src/views/Collection/hooks/useDialog";

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
  const { onClose, onOpen, open } = useDialog();

  // Fetch the dataset assets on open of download modal.
  const { datasetAssets, isError, isLoading } = useDatasetAssets(
    datasetId,
    open
  );

  return (
    <>
      <Button
        datasetName={name}
        data-testid="dataset-download-button"
        disabled={isDisabled}
        onClick={onOpen}
      />
      <StyledDialog onClose={onClose} open={open}>
        <Content
          dataAssets={datasetAssets}
          isError={isError}
          isLoading={isLoading}
          name={name}
          onClose={onClose}
        />
      </StyledDialog>
    </>
  );
};

export default DownloadDataset;
