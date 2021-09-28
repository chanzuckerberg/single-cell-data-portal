import { Classes } from "@blueprintjs/core";
import { Fragment } from "react";
import {
  Collection,
  COLLECTION_LINK_TYPE,
  COLLECTION_LINK_TYPE_OPTIONS,
  Dataset,
  DATASET_ASSET_FORMAT,
  Link,
  VISIBILITY_TYPE,
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

export function renderLinks(links: Link[]): (JSX.Element | null)[] {
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

      let urlName: string | null = name;

      if (!urlName) {
        if (linkTypeOption === COLLECTION_LINK_TYPE_OPTIONS["DOI"]) {
          urlName = new URL(url).pathname.substring(1);
        } else {
          urlName = getUrlHost(url);
        }
      }

      const { text } = linkTypeOption;

      if (!urlName) return null;

      return (
        <Fragment key={`${type}+${url}`}>
          <span className={Classes.TEXT_MUTED}>{text}</span>
          <StyledLink
            data-test-id="collection-link"
            target="_blank"
            rel="noopener"
            href={url}
          >
            {urlName}
          </StyledLink>
        </Fragment>
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
): JSX.Element | null {
  if (!contact_name || !contact_email) return null;

  return (
    <>
      <span className={Classes.TEXT_MUTED}>Contact</span>
      <StyledLink
        data-test-id="collection-contact"
        href={"mailto:" + contact_email}
      >
        {contact_name}
      </StyledLink>
    </>
  );
}

export function getIsPublishable(datasets: Array<Dataset>): boolean {
  return (
    datasets?.length > 0 &&
    datasets.every((dataset) => {
      // const numOfAssets = dataset.dataset_assets.length;
      const numOfDeployments = dataset.dataset_deployments.length;

      // TODO(seve): uncomment old check when loom is no longer served from the backend
      // return (
      //   numOfDeployments === 1 &&
      //   numOfAssets >= Object.keys(DATASET_ASSET_FORMAT).length
      // );

      const assetTypes = dataset.dataset_assets.map((asset) => asset.filetype);

      const hasAllFormats = Object.values(DATASET_ASSET_FORMAT).every(
        (format) => assetTypes.includes(format)
      );
      return numOfDeployments === 1 && hasAllFormats;
    }) &&
    datasets.some((dataset) => !dataset.tombstone)
  );
}
// check if revision meets publishing criteria
export function revisionIsPublishable(
  collection: Collection,
  revisionsEnabled: boolean
): boolean {
  if (!revisionsEnabled) return true;

  if (isPrivateRevision(collection)) {
    return collection.revision_diff;
  }

  return true;
}

function isPrivateRevision(collection: Collection) {
  return (
    collection.visibility === VISIBILITY_TYPE.PRIVATE && collection.has_revision
  );
}
