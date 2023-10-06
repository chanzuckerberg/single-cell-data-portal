import React, { FC, useEffect, useState } from "react";
import { API } from "src/common/API";
import { Dataset, DATASET_ASSET_FORMAT } from "src/common/entities";
import { DEFAULT_FETCH_OPTIONS } from "src/common/queries/common";
import { apiTemplateToUrl } from "src/common/utils/apiTemplateToUrl";
import { API_URL } from "src/configs/configs";
import CurlLink from "./components/CurlLink";
import DataFormat from "./components/DataFormat";
import Details from "./components/Details";
import Name from "./components/Name";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import {
  Button,
  DialogActions,
  DialogContent,
  DialogTitle,
} from "@czi-sds/components";
import { DialogLoader as Loader } from "src/components/Datasets/components/DownloadDataset/style";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { FEATURES } from "src/common/featureFlags/features";
import { Spinner } from "@blueprintjs/core";

interface Props {
  isError?: boolean;
  isLoading?: boolean;
  onClose: () => void;
  name: string;
  dataAssets: Dataset["dataset_assets"];
}

const Content: FC<Props> = ({
  isError = false,
  isLoading = false,
  onClose,
  name,
  dataAssets,
}) => {
  const isDownloadUX = useFeatureFlag(FEATURES.DOWNLOAD_UX);
  const [selectedFormat, setSelectedFormat] = useState<
    DATASET_ASSET_FORMAT | ""
  >("");
  const [fileSize, setFileSize] = useState<number>(0);
  const [fileName, setFileName] = useState<string>("");
  const [downloadLink, setDownloadLink] = useState<string>("");
  const [isDownloadLinkLoading, setIsDownloadLinkLoading] =
    useState<boolean>(false);
  const isDownloadDisabled =
    !downloadLink || isDownloadLinkLoading || isError || isLoading;
  const DialogLoader = isDownloadUX ? (
    <Loader sdsStyle="minimal" />
  ) : (
    <Spinner size={20} />
  ); // TODO: #5566 hidden under feature flag.

  useEffect(() => {
    if (!selectedFormat) return;

    const asset = dataAssets.filter(
      (dataAsset) => dataAsset.filetype === selectedFormat
    );

    if (!asset.length) {
      throw Error(`Format ${selectedFormat} not available`);
    }

    const { dataset_id: datasetId, filename, id: assetId } = asset[0];

    getDownloadLink({
      assetId,
      datasetId,
      filename,
      isDownloadUX,
      setFileName,
      setFileSize,
      setIsDownloadLinkLoading,
    });

    async function getDownloadLink({
      assetId,
      datasetId,
      filename,
      isDownloadUX,
      setFileName,
      setFileSize,
      setIsDownloadLinkLoading,
    }: GetDownloadLinkArgs) {
      setIsDownloadLinkLoading(true);

      const replace = {
        asset_id: assetId,
        dataset_id: datasetId,
      };

      const url = apiTemplateToUrl(API.DATASET_ASSET_DOWNLOAD_LINK, replace);

      try {
        const result = await (
          await fetch(`${API_URL}${url}`, {
            ...DEFAULT_FETCH_OPTIONS,
            method: isDownloadUX ? "GET" : "POST",
          })
        ).json();

        const { file_size } = result;
        const downloadURL = isDownloadUX ? result.url : result.presigned_url;

        setFileSize(file_size);
        setDownloadLink(downloadURL);
        setFileName(filename);
      } catch (error) {
        console.error("Please try again");
      }

      setIsDownloadLinkLoading(false);
    }
  }, [selectedFormat, dataAssets, isDownloadUX]);

  /**
   * Tracks dataset download analytics as specified by the custom analytics event.
   * @param event - Custom analytics event.
   * @param dataFormat - Data format (optionally specified with event "DOWNLOAD_DATA_FORMAT_CLICKED").
   */
  const handleAnalytics = (
    event: EVENTS,
    dataFormat?: DATASET_ASSET_FORMAT
  ) => {
    track(event, {
      dataset_name: name,
      data_format: dataFormat || selectedFormat,
    });
  };

  const handleChange = (format: DATASET_ASSET_FORMAT) => {
    setSelectedFormat(format);
    handleAnalytics(EVENTS.DOWNLOAD_DATA_FORMAT_CLICKED, format);
  };

  const availableFormats = dataAssets.map((dataAsset) => dataAsset.filetype);

  return (
    <>
      <DialogTitle title="Download Dataset" />
      <DialogContent>
        {isError && <div>Dataset download is currently not available.</div>}
        {isLoading && DialogLoader}
        {!isError && !isLoading && (
          <>
            <Name name={name} />
            <DataFormat
              availableFormats={availableFormats}
              handleChange={handleChange}
              isDisabled={isDownloadLinkLoading}
              selectedFormat={selectedFormat}
            />
            <Details
              curlPreview={
                downloadLink &&
                !isDownloadLinkLoading && (
                  <CurlLink
                    fileName={fileName}
                    handleAnalytics={() =>
                      handleAnalytics(EVENTS.DOWNLOAD_DATA_COPY)
                    }
                    link={downloadLink}
                  />
                )
              }
              fileSize={fileSize}
              isLoading={isDownloadLinkLoading}
              selected={Boolean(fileSize)}
            />
          </>
        )}
      </DialogContent>
      <DialogActions>
        <Button
          isAllCaps={false}
          onClick={onClose}
          sdsStyle="minimal"
          sdsType="secondary"
        >
          Cancel
        </Button>
        <Button
          data-testid="download-asset-download-button"
          disabled={isDownloadDisabled}
          href={downloadLink}
          onClick={() => handleAnalytics(EVENTS.DOWNLOAD_DATA_COMPLETE)}
          sdsStyle="square"
          sdsType="primary"
        >
          Download
        </Button>
      </DialogActions>
    </>
  );

  interface GetDownloadLinkArgs {
    assetId: string;
    datasetId: string;
    filename: string;
    isDownloadUX: boolean;
    setFileName: (value: string) => void;
    setFileSize: (value: number) => void;
    setIsDownloadLinkLoading: (value: boolean) => void;
  }
};

export default Content;
