import { Classes, Intent } from "@blueprintjs/core";
import { FC, useEffect, useState } from "react";
import { API } from "src/common/API";
import { Dataset, DATASET_ASSET_FORMAT } from "src/common/entities";
import { DEFAULT_FETCH_OPTIONS } from "src/common/queries/common";
import { apiTemplateToUrl } from "src/common/utils/apiTemplateToUrl";
import { API_URL } from "src/configs/configs";
import CurlLink from "./components/CurlLink";
import DataFormat from "./components/DataFormat";
import Details from "./components/Details";
import Name from "./components/Name";
import { CancelButton, DownloadButton, Wrapper } from "./style";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";

interface Props {
  onClose: () => void;
  name: string;
  dataAssets: Dataset["dataset_assets"];
}

const Content: FC<Props> = ({ onClose, name, dataAssets }) => {
  const [selectedFormat, setSelectedFormat] = useState<
    DATASET_ASSET_FORMAT | ""
  >("");
  const [fileSize, setFileSize] = useState<number>(0);
  const [fileName, setFileName] = useState<string>("");
  const [downloadLink, setDownloadLink] = useState<string>("");
  const [isLoading, setIsLoading] = useState<boolean>(false);

  useEffect(() => {
    if (!selectedFormat) return;

    const asset = dataAssets.filter(
      (dataAsset) => dataAsset.filetype === selectedFormat
    );

    if (!asset.length) {
      throw Error(`Format ${selectedFormat} not available`);
    }

    const { dataset_id: datasetId, id: assetId } = asset[0];

    getDownloadLink({
      assetId,
      datasetId,
      setFileName,
      setFileSize,
      setIsLoading,
    });

    async function getDownloadLink({
      assetId,
      datasetId,
      setFileName,
      setFileSize,
      setIsLoading,
    }: GetDownloadLinkArgs) {
      setIsLoading(true);

      const replace = {
        asset_id: assetId,
        dataset_id: datasetId,
      };

      const url = apiTemplateToUrl(API.DATASET_ASSET_DOWNLOAD_LINK, replace);

      try {
        const result = await (
          await fetch(`${API_URL}${url}`, {
            ...DEFAULT_FETCH_OPTIONS,
            method: "POST",
          })
        ).json();

        const { file_size, presigned_url, file_name } = result;

        setFileSize(file_size);
        setDownloadLink(presigned_url);
        setFileName(file_name);
      } catch (error) {
        console.error("Please try again");
      }

      setIsLoading(false);
    }
  }, [selectedFormat, dataAssets]);

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

  const renderDownload = () => {
    return (
      <DownloadButton
        disabled={!downloadLink || isLoading}
        data-test-id="download-asset-download-button"
        href={downloadLink}
        intent={Intent.PRIMARY}
        onClick={() => handleAnalytics(EVENTS.DOWNLOAD_DATA_COMPLETE)}
      >
        Download
      </DownloadButton>
    );
  };

  const availableFormats = dataAssets.map((dataAsset) => dataAsset.filetype);

  return (
    <>
      <div className={Classes.DIALOG_BODY}>
        <Wrapper>
          <Name name={name} />
          <DataFormat
            handleChange={handleChange}
            isDisabled={isLoading}
            selectedFormat={selectedFormat}
            availableFormats={availableFormats}
          />
          <Details
            isLoading={isLoading}
            fileSize={fileSize}
            selected={Boolean(fileSize)}
          />
          {downloadLink && !isLoading && (
            <CurlLink
              fileName={fileName}
              handleAnalytics={() => handleAnalytics(EVENTS.DOWNLOAD_DATA_COPY)}
              link={downloadLink}
            />
          )}
        </Wrapper>
      </div>
      <div className={Classes.DIALOG_FOOTER}>
        <CancelButton onClick={onClose} minimal>
          Cancel
        </CancelButton>
        {renderDownload()}
      </div>
    </>
  );

  interface GetDownloadLinkArgs {
    assetId: string;
    datasetId: string;
    setFileName: (value: string) => void;
    setFileSize: (value: number) => void;
    setIsLoading: (value: boolean) => void;
  }
};

export default Content;
