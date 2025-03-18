import React, { FC } from "react";
import { DATASET_ASSET_FORMAT } from "src/common/entities";
import { FormLabel } from "@mui/material";
import { FormControl } from "./style";
import { InputCheckbox } from "@czi-sds/components";
import { DownloadLinkType } from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content";
import { possibleDownloadFormats } from "./constants";
interface Props {
  handleChange: (format: DATASET_ASSET_FORMAT) => void;
  isDisabled?: boolean;
  selectedFormats: DATASET_ASSET_FORMAT[];
  availableFormats: Set<DATASET_ASSET_FORMAT>;
  downloadLinks: DownloadLinkType[];
}

const DataFormat: FC<Props> = ({
  handleChange: handleChangeRaw,
  isDisabled = false,
  selectedFormats,
  availableFormats,
  downloadLinks,
}) => {
  const handleChange = (event: React.FormEvent<HTMLElement>) => {
    const value = (event.target as HTMLInputElement)
      .value as DATASET_ASSET_FORMAT;
    handleChangeRaw(value);
  };

  const humanFileSize = (size: number) => {
    const i = size == 0 ? 0 : Math.floor(Math.log(size) / Math.log(1024));
    return (
      +(size / Math.pow(1024, i)).toFixed(0) * 1 +
      " " +
      ["B", "KB", "MB", "GB", "TB"][i]
    );
  };

  const getFileSize = (
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

  return (
    <FormControl>
      {possibleDownloadFormats.map(
        ({ format, label, type, description }) =>
          availableFormats.has(format) && (
            <div className="data-format-info" key={format}>
              <FormLabel className="data-format-type">{type}</FormLabel>
              <span className="data-format-checkbox-group">
                <InputCheckbox
                  disabled={isDisabled || !availableFormats.has(format)}
                  label={
                    <span>
                      {label}
                      <span className="file-size">
                        {getFileSize(downloadLinks, format)}
                      </span>
                    </span>
                  }
                  value={format}
                  onChange={handleChange}
                  checked={selectedFormats.includes(format)}
                />
              </span>
              <p className="data-type-description">{description}</p>
            </div>
          )
      )}
    </FormControl>
  );
};

export default DataFormat;
