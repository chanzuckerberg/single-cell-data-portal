/* eslint-disable @typescript-eslint/camelcase */
import React, { FC } from "react";
import { DatasetDeployment } from "src/common/entities";
import { StyledAnchor, Wrapper } from "./style";

interface Props {
  deployments: DatasetDeployment[];
}

const ViewDataset: FC<Props> = ({ deployments }) => {
  const deployment = deployments[0];

  return (
    <Wrapper>
      <StyledAnchor
        key={deployment.link}
        href={deployment.link}
        target="_blank"
        rel="noopener"
      >
        Remix
        <br />
      </StyledAnchor>
    </Wrapper>
  );
};

export default ViewDataset;
