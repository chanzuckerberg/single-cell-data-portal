import { DATASET_ASSET_FORMAT } from "src/common/entities";
import { DownloadLinkType } from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content";

export const humanFileSize = (size: number) => {
  const i = size == 0 ? 0 : Math.floor(Math.log(size) / Math.log(1024));
  return (
    +(size / Math.pow(1024, i)).toFixed(0) * 1 +
    " " +
    ["B", "KB", "MB", "GB", "TB"][i]
  );
};

export const getFileSize = (
  downloadLinks: DownloadLinkType[],
  fileType: DATASET_ASSET_FORMAT
) => {
  let fileSize: number | string = "--";
  if (downloadLinks.length === 0) {
    return fileSize;
  }
  if (fileType === DATASET_ASSET_FORMAT.ATAC_INDEX) {
    fileSize = downloadLinks
      .filter(
        (x) =>
          x.filetype === DATASET_ASSET_FORMAT.ATAC_INDEX ||
          x.filetype === DATASET_ASSET_FORMAT.ATAC_FRAGMENT
      )
      .reduce((acc, x) => {
        return x.fileSize ? acc + x.fileSize : acc;
      }, 0);
  } else {
    fileSize =
      downloadLinks.find((x) => x.filetype === fileType)?.fileSize || "--";
  }
  if (typeof fileSize === "string") {
    return fileSize;
  }
  return humanFileSize(fileSize);
};
