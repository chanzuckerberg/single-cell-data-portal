import { Classes } from "@blueprintjs/core";
import React from "react";
import {
  Collection,
  COLLECTION_LINK_TYPE,
  COLLECTION_LINK_TYPE_OPTIONS,
  DATASET_ASSET_FORMAT,
  Link,
} from "src/common/entities";
import { getUrlHost } from "src/common/utils/getUrlHost";
import { StyledLink } from "./common/style";

const LINK_ORDER: COLLECTION_LINK_TYPE[] = [
  COLLECTION_LINK_TYPE.DOI,
  COLLECTION_LINK_TYPE.DATA_SOURCE,
  COLLECTION_LINK_TYPE.RAW_DATA,
  COLLECTION_LINK_TYPE.PROTOCOL,
  COLLECTION_LINK_TYPE.LAB_WEBSITE,
  COLLECTION_LINK_TYPE.OTHER,
];

export function renderLinks(links: Link[]) {
  const linkTypesToLinks = createLinkTypesToLinks(links);

  const orderedLinks: Link[] = [];

  for (const linkType of LINK_ORDER) {
    const links = linkTypesToLinks[linkType];

    if (!links) continue;

    for (const link of links) {
      orderedLinks.push(link);
    }
  }

  return orderedLinks.map(
    ({ link_url: url, link_type: type, link_name: name }) => {
      const linkTypeOption = COLLECTION_LINK_TYPE_OPTIONS[type];

      if (!linkTypeOption) return null;

      const urlName = name ? name : getUrlHost(url);

      const { text } = linkTypeOption;

      if (!urlName) return null;

      return (
        <React.Fragment key={`${type}+${url}`}>
          <span className={Classes.TEXT_MUTED}>{text}</span>
          <StyledLink href={url}>{urlName}</StyledLink>
        </React.Fragment>
      );
    }
  );
}

type LinkTypesToLinks = {
  [key in COLLECTION_LINK_TYPE]?: Link[];
};

function createLinkTypesToLinks(links: Link[]) {
  const linkTypesToLinks: LinkTypesToLinks = {};

  for (const link of links) {
    const linkTypeToLinks = linkTypesToLinks[link.link_type] || [];

    linkTypeToLinks.push(link);

    linkTypesToLinks[link.link_type] = linkTypeToLinks;
  }

  return linkTypesToLinks;
}

export function renderContact(
  contact_name?: Collection["contact_name"],
  contact_email?: Collection["contact_email"]
) {
  if (!contact_name || !contact_email) return null;

  return (
    <>
      <span className={Classes.TEXT_MUTED}>Contact</span>
      <StyledLink href={"mailto:" + contact_email}>{contact_name}</StyledLink>
    </>
  );
}

export function getIsPublishable(datasets: Collection["datasets"]) {
  return (
    Boolean(datasets?.length) &&
    datasets.every((dataset) => {
      const numOfAssets = dataset.dataset_assets.length;
      const numOfDeployments = dataset.dataset_deployments.length;

      return (
        numOfAssets + numOfDeployments ===
        Object.keys(DATASET_ASSET_FORMAT).length
      );
    })
  );
}
