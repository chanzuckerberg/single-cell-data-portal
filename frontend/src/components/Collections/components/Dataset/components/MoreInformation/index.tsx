import React, { FC } from "react";
import { Collection } from "src/common/entities";
import { SmallColumn } from "../../common/style";
import { StyledAnchor, Wrapper } from "./style";
interface AnchorProps {
  url: string;
}

const Anchor: FC<AnchorProps> = ({ url, children }) => {
  return (
    <StyledAnchor href={url} target="_blank" rel="noopener">
      {children}
    </StyledAnchor>
  );
};

interface Props {
  links: Collection["links"];
}

const MoreInformation: FC<Props> = ({ links }) => {
  const uniqueLinks = Array.from(
    new Map(links.map((link) => [link.url, link.name]))
  );

  const styledLinks = uniqueLinks.map((link, index) => {
    const [url, name] = link;

    if (index === uniqueLinks.length - 1) {
      return (
        <div key={url}>
          <Anchor url={url}>{name}</Anchor>
        </div>
      );
    } else {
      return (
        <div key={url}>
          <Anchor url={url}>{name},&nbsp;</Anchor>
        </div>
      );
    }
  });

  return (
    <SmallColumn>
      <Wrapper>{styledLinks}</Wrapper>
    </SmallColumn>
  );
};

export default MoreInformation;
