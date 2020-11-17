import React, { FC } from "react";
import { Collection } from "src/common/entities";
import { getUrlHost } from "src/common/utils/getUrlHost";
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
    // eslint-disable-next-line @typescript-eslint/camelcase
    new Map(links.map(({ link_url, link_name }) => [link_url, link_name]))
  );

  const styledLinks = uniqueLinks.map((link, index) => {
    const [url, name] = link;

    const displayName = name || getUrlHost(url);

    if (index === uniqueLinks.length - 1) {
      return (
        <div key={url}>
          <Anchor url={url}>{displayName}</Anchor>
        </div>
      );
    } else {
      return (
        <div key={url}>
          <Anchor url={url}>{displayName},&nbsp;</Anchor>
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
