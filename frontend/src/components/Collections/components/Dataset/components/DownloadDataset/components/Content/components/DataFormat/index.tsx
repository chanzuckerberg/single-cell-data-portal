import React, { FC } from "react";
import { DATASET_ASSET_FORMAT } from "src/common/entities";
import { FormControl, FormLabel, RadioGroup } from "@mui/material";
import { InputRadio, Tooltip } from "@czi-sds/components";
import {
  DOWNLOAD_TOOLTIP_SLOT_PROPS,
  DOWNLOAD_TOOLTIP_TITLE,
  TOOLTIP_SLOT_PROPS,
  TOOLTIP_TITLE,
} from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/DataFormat/constants";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { FEATURES } from "src/common/featureFlags/features";

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
  const isDownloadUX = useFeatureFlag(FEATURES.DOWNLOAD_UX);
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
          slotProps={
            isDownloadUX ? TOOLTIP_SLOT_PROPS : DOWNLOAD_TOOLTIP_SLOT_PROPS
          } // TODO(cc) Download UI #5566 hidden under feature flag.
          title={
            isRDS ? null : isDownloadUX ? TOOLTIP_TITLE : DOWNLOAD_TOOLTIP_TITLE
          } // TODO(cc) Download UI #5566 hidden under feature flag.
        >
          {/* The radio button is enclosed within a <span> tag to enable tooltip functionality when the radio button is disabled. */}
          {/* See https://github.com/chanzuckerberg/sci-components/blob/main/packages/components/src/core/Tooltip/index.tsx#L28. */}
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
