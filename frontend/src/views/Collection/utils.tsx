import { Classes } from "@blueprintjs/core";
import React from "react";
import {
  Collection,
  COLLECTION_LINK_TYPE_OPTIONS,
  DATASET_ASSET_FORMAT,
  Link,
} from "src/common/entities";
import { getUrlHost } from "src/common/utils/getUrlHost";
import { StyledLink } from "./common/style";

export function renderLinks(links: Link[]) {
  return links?.map(({ link_url: url, link_type: type }) => {
    const linkTypeOption = COLLECTION_LINK_TYPE_OPTIONS[type];

    if (!linkTypeOption) return null;

    const urlHost = getUrlHost(url);

    const { text } = linkTypeOption;

    if (!urlHost) return null;

    return (
      <React.Fragment key={`${type}+${url}`}>
        <span className={Classes.TEXT_MUTED}>{text}</span>
        <StyledLink href={url}>{urlHost}</StyledLink>
      </React.Fragment>
    );
  });
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
