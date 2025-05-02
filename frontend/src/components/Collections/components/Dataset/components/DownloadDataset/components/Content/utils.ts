import { DATASET_ASSET_FORMAT } from "src/common/entities";
import { DownloadLinkType } from "./index";
import { DEFAULT_FETCH_OPTIONS } from "src/common/queries/common";
import { apiTemplateToUrl } from "src/common/utils/apiTemplateToUrl";
import { API_URL } from "src/configs/configs";
import { API } from "src/common/API";
import { DatasetAsset } from "src/common/entities";
import { getNotAvailableText } from "./components/BrowserTab/components/DataFormat/constants";

export const downloadMultipleFiles = async (
  formatsToDownload: DATASET_ASSET_FORMAT[],
  downloadLinks: DownloadLinkType[]
) => {
  for (const format of formatsToDownload) {
    const downloadLink = downloadLinks.find((link) => link.filetype === format);

    if (
      downloadLink &&
      downloadLink.downloadURL &&
      downloadLink.downloadURL !== getNotAvailableText(downloadLink.filetype)
    ) {
      const a = document.createElement("a");
      a.href = downloadLink.downloadURL;
      a.download = downloadLink.downloadURL.split("/").pop() || "";
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);

      // Wait for 500ms before downloading the next file
      // This is to prevent the browser from blocking the downloads
      await new Promise((resolve) => setTimeout(resolve, 500));
    }
  }
};

export const getDownloadLink = async (
  asset: DatasetAsset
): Promise<DownloadLinkType | null> => {
  try {
    const url = apiTemplateToUrl(API.DATASET_ASSET_DOWNLOAD_LINK, {
      asset_id: asset.id,
      dataset_id: asset.dataset_id,
    });

    const result = await fetch(`${API_URL}${url}`, {
      ...DEFAULT_FETCH_OPTIONS,
      method: "GET",
    }).then((res) => res.json());
    return {
      filetype: asset.filetype,
      fileSize: result.file_size,
      downloadURL: result.url,
    };
  } catch (error) {
    console.error("Error fetching download link for asset:", asset, error);
    return null;
  }
};

export enum ActiveTab {
  Browser = "Browser",
  PythonApi = "PythonApi",
  RApi = "RApi",
}
