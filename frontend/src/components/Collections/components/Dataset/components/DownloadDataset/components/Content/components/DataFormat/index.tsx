import {
  Intent,
  PopoverInteractionKind,
  Radio,
  RadioGroup,
  Tooltip,
} from "@blueprintjs/core";
import * as React from "react";
import { FC } from "react";
import { DATASET_ASSET_FORMAT } from "src/common/entities";
import { Section, Title } from "../common/style";

interface Props {
  handleChange: (format: DATASET_ASSET_FORMAT) => void;
  isDisabled: boolean;
  selectedFormat: DATASET_ASSET_FORMAT | "";
  availableFormats: DATASET_ASSET_FORMAT[];
  rdsSkipped: boolean;
}

const DataFormat: FC<Props> = ({
  handleChange: handleChangeRaw,
  isDisabled = false,
  selectedFormat,
  availableFormats,
  rdsSkipped,
}) => {
  const handleChange = (event: React.FormEvent<HTMLElement>) => {
    const value = (event.target as HTMLInputElement)
      .value as DATASET_ASSET_FORMAT;

    handleChangeRaw(value);
  };

  return (
    <Section>
      <Title>DATA FORMAT</Title>
      <RadioGroup
        inline
        name="dataFormat"
        disabled={isDisabled}
        onChange={handleChange}
        selectedValue={selectedFormat}
      >
        <Radio
          disabled={!availableFormats.includes(DATASET_ASSET_FORMAT.H5AD)}
          label=".h5ad (AnnData v0.7)"
          value={DATASET_ASSET_FORMAT.H5AD}
        />
        <Tooltip
          disabled={!rdsSkipped}
          interactionKind={PopoverInteractionKind.HOVER}
          content="A .rds (Seurat v3) download is unavailable due to limitations in the R dgCMatrix sparse matrix class."
          intent={Intent.DANGER}
        >
          <Radio
            disabled={!availableFormats.includes(DATASET_ASSET_FORMAT.RDS)}
            label=".rds (Seurat v3)"
            value={DATASET_ASSET_FORMAT.RDS}
          />
        </Tooltip>
      </RadioGroup>
    </Section>
  );
};

export default DataFormat;
