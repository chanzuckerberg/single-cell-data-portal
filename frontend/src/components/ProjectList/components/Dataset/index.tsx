import React, { FC } from "react";
import { Dataset as IDataset, Project as IProject } from "src/common/entities";
import MoreInformation from "./components/MoreInformation";
import ViewDataset from "./components/ViewDataset";
import { Name, Wrapper } from "./style";

interface Props {
  dataset: IDataset;
  links: IProject["links"];
}

const Dataset: FC<Props> = ({ dataset, links }) => {
  return (
    <Wrapper>
      <Name>{dataset.name}</Name>
      <ViewDataset deployments={dataset.dataset_deployments} />
      <MoreInformation links={links} />
    </Wrapper>
  );
};

export default Dataset;
