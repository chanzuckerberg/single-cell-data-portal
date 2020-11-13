/* eslint-disable @typescript-eslint/camelcase */
import { Classes } from "@blueprintjs/core";
import React, { FC, useEffect, useState } from "react";
import { API } from "src/common/API";
import { Dataset, DATASET_ASSET_FORMAT } from "src/common/entities";
import { apiTemplateToUrl } from "src/common/utils/apiTemplateToUrl";
import { API_URL } from "src/configs/configs";
import CurlLink from "./components/CurlLink";
import DataFormat from "./components/DataFormat";
import Details from "./components/Details";
import Name from "./components/Name";
import { Cancel, DisabledDownload, Download, Wrapper } from "./style";

interface Props {
  onClose: () => void;
  name: string;
  dataAssets: Dataset["dataset_assets"];
}

const Content: FC<Props> = ({ onClose, name, dataAssets }) => {
  const [format, setFormat] = useState<DATASET_ASSET_FORMAT | "">("");
  const [fileSize, setFileSize] = useState<number>(0);
  const [fileName, setFileName] = useState<string>("");
  const [downloadLink, setDownloadLink] = useState<string | null>(null);
  const [isLoading, setIsLoading] = useState<boolean>(false);

  useEffect(() => {
    if (!format) return;

    const asset = dataAssets.filter(
      (dataAsset) => dataAsset.filetype === format
    );

    if (!asset.length) {
      throw Error(`Format ${format} not available`);
    }

    const { dataset_id: datasetId, id: assetId } = asset[0];

    getDownloadLink({
      assetId,
      datasetId,
      setFileName,
      setFileSize,
      setIsLoading,
    });
  }, [format, dataAssets]);

  const handleChange = (format: DATASET_ASSET_FORMAT) => {
    setFormat(format);
  };

  const renderDownload = () => {
    if (!downloadLink || isLoading) {
      return <DisabledDownload>Download</DisabledDownload>;
    }

    return (
      <Download
        data-test-id="download-asset-download-button"
        href={downloadLink}
      >
        Download
      </Download>
    );
  };

  return (
    <>
      <div className={Classes.DIALOG_BODY}>
        <Wrapper>
          <Name name={name} />
          <DataFormat
            handleChange={handleChange}
            isDisabled={isLoading}
            format={format}
          />
          <Details
            isLoading={isLoading}
            fileSize={fileSize}
            selected={Boolean(fileSize)}
          />
          {downloadLink && !isLoading && (
            <CurlLink fileName={fileName} link={downloadLink} />
          )}
        </Wrapper>
      </div>
      <div className={Classes.DIALOG_FOOTER}>
        <Cancel onClick={onClose}>Cancel</Cancel>
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

  async function getDownloadLink({
    assetId,
    datasetId,
    setFileName,
    setFileSize,
    setIsLoading,
  }: GetDownloadLinkArgs) {
    setIsLoading(true);

    const replace = {
      asset_uuid: assetId,
      dataset_uuid: datasetId,
    };

    const url = apiTemplateToUrl(API.DATASET_ASSET_DOWNLOAD_LINK, replace);

    try {
      const result = await (
        await fetch(`${API_URL}${url}`, { method: "POST" })
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
};

export default Content;
