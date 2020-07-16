import React, { FC } from "react";
import { Project } from "src/common/entities";
import { StyledAnchor, Wrapper } from "./style";

interface Props {
  links: Project["links"];
}

const MoreInformation: FC<Props> = ({ links }) => {
  return (
    <Wrapper>
      {links.map(link => {
        return (
          <StyledAnchor
            key={link.url}
            href={link.url}
            target="_blank"
            rel="noopener"
          >
            {link.name}
            <br />
          </StyledAnchor>
        );
      })}
    </Wrapper>
  );
};

export default MoreInformation;
