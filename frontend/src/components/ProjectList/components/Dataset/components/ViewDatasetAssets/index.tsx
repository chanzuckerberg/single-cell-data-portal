/* eslint-disable @typescript-eslint/camelcase */
import React, { ChangeEvent, FC, useEffect, useState } from "react";
import { DatasetAsset as IDatasetAsset } from "src/common/entities";
import DatasetAsset from "./DatasetAsset";
import { Arrow, SelectWrapper, StyledSelect, Wrapper } from "./style";

const OPTION_INSTRUCTION = "Explore in cellxgene";

interface Props {
  datasetAssets: IDatasetAsset[];
}

const ViewDatasetAssets: FC<Props> = ({ datasetAssets }) => {
  const [link, setLink] = useState<string>("");

  useEffect(() => {
    if (!link || link === OPTION_INSTRUCTION) return;

    window.open(link, "_blank");
  }, [link]);

  const handleChange = (event: ChangeEvent<HTMLSelectElement>) =>
    setLink(event.target.value);

  return (
    <Wrapper>
      <SelectWrapper>
        <StyledSelect onChange={handleChange}>
          <option>{OPTION_INSTRUCTION}</option>
          {datasetAssets.map(datasetAsset => {
            return (
              <DatasetAsset key={datasetAsset.id} datasetAsset={datasetAsset} />
            );
          })}
        </StyledSelect>
        <Arrow />
      </SelectWrapper>
    </Wrapper>
  );
};

export default ViewDatasetAssets;
