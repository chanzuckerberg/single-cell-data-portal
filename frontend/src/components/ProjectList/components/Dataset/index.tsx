import React, { FC } from "react";
import { Dataset as IDataset, Project as IProject } from "src/common/entities";
import MoreInformation from "./components/MoreInformation";
import ViewDatasetAssets from "./components/ViewDatasetAssets";
import { Name, Wrapper } from "./style";

interface Props {
  dataset: IDataset;
  links: IProject["links"];
}

const Dataset: FC<Props> = ({ dataset, links }) => {
  return (
    <Wrapper>
      <Name>{dataset.name}</Name>
      <ViewDatasetAssets datasetAssets={dataset.dataset_assets} />
      <MoreInformation links={links} />
    </Wrapper>
  );
};

export default Dataset;
