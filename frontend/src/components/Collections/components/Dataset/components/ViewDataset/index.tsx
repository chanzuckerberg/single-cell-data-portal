/* eslint-disable @typescript-eslint/camelcase */
import React, { FC } from "react";
import { DatasetDeployment } from "src/common/entities";
import { SmallColumn } from "../../common/style";
import { StyledAnchor } from "./style";

interface Props {
  deployments: DatasetDeployment[];
}

const ViewDataset: FC<Props> = ({ deployments }) => {
  // (thuang): Temp
  // Currently BE only returns 0 or 1 Remix cellxgene link
  const deployment = deployments[0];

  return <SmallColumn>{deployment ? <View /> : "N/A"}</SmallColumn>;

  function View() {
    return (
      <StyledAnchor
        key={deployment.url}
        href={deployment.url}
        target="_blank"
        rel="noopener"
        data-test-id="view-dataset-link"
      >
        View
        <br />
      </StyledAnchor>
    );
  }
};

export default ViewDataset;
