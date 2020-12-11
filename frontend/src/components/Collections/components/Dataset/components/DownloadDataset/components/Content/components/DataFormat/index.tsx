import { Radio, RadioGroup } from "@blueprintjs/core";
import React, { FC } from "react";
import { DATASET_ASSET_FORMAT } from "src/common/entities";
import { Section, Title } from "../common/style";

interface Props {
  handleChange: (format: DATASET_ASSET_FORMAT) => void;
  isDisabled: boolean;
  format: DATASET_ASSET_FORMAT | "";
}

const DataFormat: FC<Props> = ({
  handleChange: handleChangeRaw,
  isDisabled = false,
  format,
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
        selectedValue={format}
      >
        <Radio label=".h5ad (AnnData v0.7)" value={DATASET_ASSET_FORMAT.H5AD} />
        <Radio label=".loom" value={DATASET_ASSET_FORMAT.LOOM} />
        <Radio label=".rds (Seurat v3)" value={DATASET_ASSET_FORMAT.RDS} />
      </RadioGroup>
    </Section>
  );
};

export default DataFormat;
