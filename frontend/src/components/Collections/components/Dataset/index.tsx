import React, { FC } from "react";
import {
  Collection as ICollection,
  Dataset as IDataset,
} from "src/common/entities";
import { hasAssets } from "src/common/modules/datasets/selectors";
import { SmallColumn } from "./common/style";
import DownloadDataset from "./components/DownloadDataset";
import MoreInformation from "./components/MoreInformation";
import ViewDataset from "./components/ViewDataset";
import { Name, NameChild, Wrapper } from "./style";

interface Props {
  dataset: IDataset;
  links: ICollection["links"];
}

const Dataset: FC<Props> = ({ dataset, links }) => {
  return (
    <Wrapper>
      <Name data-test-id="dataset-name" title={dataset.name}>
        <NameChild>{dataset.name}</NameChild>
      </Name>
      <ViewDataset deployments={dataset.dataset_deployments} />

      <SmallColumn>
        <DownloadDataset
          isDisabled={!hasAssets(dataset)}
          name={dataset.name}
          dataAssets={dataset.dataset_assets}
        />
      </SmallColumn>

      <MoreInformation links={links} />
    </Wrapper>
  );
};

export default Dataset;
