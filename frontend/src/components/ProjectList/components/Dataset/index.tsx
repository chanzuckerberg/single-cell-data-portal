import React, { FC } from "react";
import { Dataset as IDataset, Project as IProject } from "src/common/entities";
import DownloadDataset from "./components/DownloadDataset";
import MoreInformation from "./components/MoreInformation";
import ViewDataset from "./components/ViewDataset";
import { Name, NameChild, Wrapper } from "./style";

interface Props {
  dataset: IDataset;
  links: IProject["links"];
}

const Dataset: FC<Props> = ({ dataset, links }) => {
  return (
    <Wrapper>
      <Name data-test-id="dataset-name" title={dataset.name}>
        <NameChild>{dataset.name}</NameChild>
      </Name>
      <ViewDataset deployments={dataset.dataset_deployments} />
      <DownloadDataset
        name={dataset.name}
        dataAssets={dataset.dataset_assets}
      />
      <MoreInformation links={links} />
    </Wrapper>
  );
};

export default Dataset;
