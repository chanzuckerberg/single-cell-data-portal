import React, { FC } from "react";
import { DATASET_ASSET_FORMAT } from "src/common/entities";
import { FormLabel } from "@mui/material";
import { FormControl } from "./style";
import { InputCheckbox } from "@czi-sds/components";
import { DownloadLinkType } from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content";
import { POSSIBLE_DOWNLOAD_FORMATS } from "./constants";
import { getFileSize } from "./utils";
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
  return (
    <FormControl>
      {POSSIBLE_DOWNLOAD_FORMATS.map(
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
