import { Spinner } from "@blueprintjs/core";
import * as React from "react";
import { FC, useCallback } from "react";
import { useDatasetAssets } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DownloadDataset/util";
import Content from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content";
import { StyledButton } from "src/components/Collections/components/Dataset/components/DownloadDataset/style";
import Modal from "src/components/common/Modal";
import { ModalContentWrapper } from "./style";

interface Props {
  Button?: React.ElementType;
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
  Button = StyledButton,
  datasetId,
  isDisabled = false,
  name,
}) => {
  // Open state of download modal.
  const [isOpen, setIsOpen] = React.useState(false);

  // Function toggling open state of modal.
  const toggleOpen = useCallback(() => {
    setIsOpen(!isOpen);
  }, [isOpen]);

  // Fetch the dataset assets on open of download modal.
  const { datasetAssets, isError, isLoading } = useDatasetAssets(
    datasetId,
    isOpen
  );

  return (
    <>
      <Button
        datasetName={name}
        disabled={isDisabled}
        onClick={toggleOpen}
        data-testid="dataset-download-button"
      />
      <Modal
        title="Download Dataset"
        isCloseButtonShown={false}
        isOpen={isOpen}
        onClose={toggleOpen}
        className={isLoading || isError ? "modal-loading" : undefined}
      >
        {isLoading ? (
          <ModalContentWrapper>
            <Spinner size={20} />
          </ModalContentWrapper>
        ) : null}
        {isError ? (
          <ModalContentWrapper>
            Dataset download is currently not available.
          </ModalContentWrapper>
        ) : null}
        {!isLoading && !isError ? (
          <Content
            name={name}
            dataAssets={datasetAssets}
            onClose={toggleOpen}
          />
        ) : null}
      </Modal>
    </>
  );
};

export default DownloadDataset;
