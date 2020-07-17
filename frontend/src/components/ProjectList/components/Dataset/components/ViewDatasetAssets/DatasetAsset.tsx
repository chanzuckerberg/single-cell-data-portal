/* eslint-disable @typescript-eslint/camelcase */
import React, { FC } from "react";
import {
  DatasetAsset as IDatasetAsset,
  DATASET_ASSET_TYPE,
} from "src/common/entities";

interface Props {
  datasetAsset: IDatasetAsset;
}

const DATASET_ASSET_TEXT = {
  [DATASET_ASSET_TYPE.ORIGINAL]: "Original dataset",
  [DATASET_ASSET_TYPE.REMIX]: "Corpora-transformed dataset",
};

const DOMAIN_PREFIX = "https://cellxgene.cziscience.com/d";

const DatasetAsset: FC<Props> = ({ datasetAsset }) => {
  const { type, format, file_name } = datasetAsset;

  const text = DATASET_ASSET_TEXT[type];

  const value = `${DOMAIN_PREFIX}/${file_name}.${format}`;

  return <option value={value}>{text}</option>;
};

export default DatasetAsset;
