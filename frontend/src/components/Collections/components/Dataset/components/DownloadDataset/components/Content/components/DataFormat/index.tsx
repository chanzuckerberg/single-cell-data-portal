import React, { FC } from "react";
import { DATASET_ASSET_FORMAT } from "src/common/entities";
import { FormControl, FormLabel, RadioGroup } from "@mui/material";
import { InputRadio, Tooltip } from "@czi-sds/components";
import {
  TOOLTIP_SLOT_PROPS,
  TOOLTIP_TITLE,
} from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/DataFormat/constants";

interface Props {
  handleChange: (format: DATASET_ASSET_FORMAT) => void;
  isDisabled: boolean;
  selectedFormat: DATASET_ASSET_FORMAT | "";
  availableFormats: DATASET_ASSET_FORMAT[];
}

const DataFormat: FC<Props> = ({
  handleChange: handleChangeRaw,
  isDisabled = false,
  selectedFormat,
  availableFormats,
}) => {
  const isH5AD = availableFormats.includes(DATASET_ASSET_FORMAT.H5AD);
  const isRDS = availableFormats.includes(DATASET_ASSET_FORMAT.RDS);
  const handleChange = (event: React.FormEvent<HTMLElement>) => {
    const value = (event.target as HTMLInputElement)
      .value as DATASET_ASSET_FORMAT;

    handleChangeRaw(value);
  };

  return (
    <FormControl>
      <FormLabel>Data Format</FormLabel>
      <RadioGroup
        name="dataFormat"
        onChange={handleChange}
        value={selectedFormat}
      >
        <InputRadio
          disabled={isDisabled || !isH5AD}
          label=".h5ad (AnnData v0.8)"
          value={DATASET_ASSET_FORMAT.H5AD}
        />
        <Tooltip
          arrow
          placement="top"
          sdsStyle="dark"
          slotProps={TOOLTIP_SLOT_PROPS}
          title={isRDS ? null : TOOLTIP_TITLE}
        >
          <span>
            <InputRadio
              disabled={isDisabled || !isRDS}
              label=".rds (Seurat v4)"
              value={DATASET_ASSET_FORMAT.RDS}
            />
          </span>
        </Tooltip>
      </RadioGroup>
    </FormControl>
  );
};

export default DataFormat;
