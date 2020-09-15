import React, { FC } from "react";
import { Project } from "src/common/entities";
import { SmallColumn } from "../../common/style";
import { StyledAnchor } from "./style";
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
  links: Project["links"];
}

const MoreInformation: FC<Props> = ({ links }) => {
  const uniqueLinks = Array.from(
    new Map(links.map((link) => [link.url, link.name]))
  );

  const styledLinks = uniqueLinks.map((link, index) => {
    const [url, name] = link;

    if (index === uniqueLinks.length - 1) {
      return (
        <Anchor key={url} url={url}>
          {name}
        </Anchor>
      );
    } else {
      return (
        <Anchor key={url} url={url}>
          {name},&nbsp;
        </Anchor>
      );
    }
  });

  return <SmallColumn>{styledLinks}</SmallColumn>;
};

export default MoreInformation;
