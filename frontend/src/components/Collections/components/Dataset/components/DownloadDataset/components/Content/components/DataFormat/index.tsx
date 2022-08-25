import {
  Intent,
  PopoverInteractionKind,
  Position,
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
}

const DataFormat: FC<Props> = ({
  handleChange: handleChangeRaw,
  isDisabled = false,
  selectedFormat,
  availableFormats,
}) => {
  const handleChange = (event: React.FormEvent<HTMLElement>) => {
    const value = (event.target as HTMLInputElement)
      .value as DATASET_ASSET_FORMAT;

    handleChangeRaw(value);
  };

  const renderH5adRadio = (): React.ReactElement => {
    return (
      <Radio
        disabled={!availableFormats.includes(DATASET_ASSET_FORMAT.H5AD)}
        label=".h5ad (AnnData v0.8)"
        value={DATASET_ASSET_FORMAT.H5AD}
      />
    );
  };

  const renderRdsRadio = (): React.ReactElement => {
    return (
      <Radio
        disabled={!availableFormats.includes(DATASET_ASSET_FORMAT.RDS)}
        label=".rds (Seurat v3)"
        value={DATASET_ASSET_FORMAT.RDS}
      />
    );
  };

  const renderDisabledRdsRadio = (): React.ReactElement => {
    return (
      <Tooltip
        disabled={false}
        interactionKind={PopoverInteractionKind.HOVER}
        content="A .rds (Seurat v3) download is unavailable due to limitations in the R dgCMatrix sparse matrix class."
        intent={Intent.DANGER}
        position={Position.TOP}
      >
        {/* Logically renders only when rds format not available, but a Tooltip wrapper also happens to prevent 
        proper functioning of radio buttons with a larger radio group outside of the Tooltip wrapper */}
        {renderRdsRadio()}
      </Tooltip>
    );
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
        {renderH5adRadio()}
        {availableFormats.includes(DATASET_ASSET_FORMAT.RDS)
          ? renderRdsRadio()
          : renderDisabledRdsRadio()}
      </RadioGroup>
    </Section>
  );
};

export default DataFormat;
