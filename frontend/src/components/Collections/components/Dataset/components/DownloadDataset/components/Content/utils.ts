import { DATASET_ASSET_FORMAT } from "src/common/entities";
import { DownloadLinkType } from "./index";
import { DEFAULT_FETCH_OPTIONS } from "src/common/queries/common";
import { apiTemplateToUrl } from "src/common/utils/apiTemplateToUrl";
import { API_URL } from "src/configs/configs";
import { API } from "src/common/API";
import { DatasetAsset } from "src/common/entities";

export const downloadMultipleFiles = (
  formatsToDownload: DATASET_ASSET_FORMAT[],
  downloadLinks: DownloadLinkType[]
) => {
  downloadLinks.forEach((download) => {
    if (formatsToDownload.includes(download.filetype) && download.downloadURL) {
      const a = document.createElement("a");
      a.href = download.downloadURL;
      a.download = download.downloadURL.split("/").pop() || "";
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
    }
  });
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
