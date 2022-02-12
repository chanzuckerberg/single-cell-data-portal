import { Classes } from "@blueprintjs/core";
import { Fragment } from "react";
import {
  Collection,
  COLLECTION_LINK_TYPE,
  COLLECTION_LINK_TYPE_OPTIONS,
  Dataset,
  DATASET_ASSET_FORMAT,
  Link,
  PROCESSING_STATUS,
  VISIBILITY_TYPE,
} from "src/common/entities";
import { getUrlHost } from "src/common/utils/getUrlHost";
import { CollectionMetadataLink } from "src/components/Collection/components/CollectionMetadata";
import { StyledLink } from "./common/style";

const LINK_ORDER: COLLECTION_LINK_TYPE[] = [
  COLLECTION_LINK_TYPE.DOI,
  COLLECTION_LINK_TYPE.DATA_SOURCE,
  COLLECTION_LINK_TYPE.RAW_DATA,
  COLLECTION_LINK_TYPE.PROTOCOL,
  COLLECTION_LINK_TYPE.LAB_WEBSITE,
  COLLECTION_LINK_TYPE.OTHER,
];

/**
 * Returns collection metadata in preferred order of display.
 * @param links - links associated with collection.
 * @param contactName - Name of collection contact.
 * @param contactEmail - Email of collection contact.
 * @param summaryCitation - Summary citation format of collection publication metadata.
 * @returns Array of collection metadata in preferred order of display.
 */
export function buildCollectionMetadataLinks(
  links: Link[],
  contactName?: Collection["contact_name"],
  contactEmail?: Collection["contact_email"],
  summaryCitation?: string
): CollectionMetadataLink[] {
  const collectionMetadataLinks = [];

  /* Clone the links so we can safely remove the DOI link type if present and display it separately from the other link
    types at the top of the metadata links list.*/
  const linksClone = [...links];

  /* If collection has an associated DOI, display either the summary citation or the DOI itself. */
  const doiLink = links.find((link: Link) => link.link_type === "DOI");
  if (doiLink) {
    const doiMetadataLink = buildDoiMetadataLink(doiLink, summaryCitation);
    if (doiMetadataLink) {
      collectionMetadataLinks.push(doiMetadataLink);

      /* Remove the DOI from links so we don't display it again in the links section below contact. */
      linksClone.splice(linksClone.indexOf(doiLink), 1);
    }
  }

  /* Add contact name and email to the top of collection metadata list. */
  if (contactName && contactEmail) {
    collectionMetadataLinks.push({
      label: "Contact",
      testId: "collection-contact",
      url: `mailto:${contactEmail}`,
      value: contactName,
    });
  }

  /* Add collection links - except DOI - to collection metadata list. */
  /* Collection metadata reordered for preferred order of display. */
  const orderedLinks = sortCollectionLinks(linksClone);
  for (const orderedLink of orderedLinks) {
    /* Build and add any valid collection metadata. */
    const collectionMetadataLink = buildCollectionMetadataLink(orderedLink);
    if (!collectionMetadataLink) continue;
    collectionMetadataLinks.push(collectionMetadataLink);
  }

  /* Return valid collection metadata. */
  return collectionMetadataLinks;
}

/**
 * @deprecated Remove once feature flag is removed (#1718).
 * Returns collection metadata in preferred order of display.
 * @param links - links associated with collection.
 * @param contactName
 * @param contactEmail
 * @returns Array of collection metadata in preferred order of display.
 */
export function buildCollectionMetadataLinksDeprecated(
  links: Link[],
  contactName?: Collection["contact_name"],
  contactEmail?: Collection["contact_email"]
): CollectionMetadataLink[] {
  const collectionMetadataLinks = [];

  /* Add contact name and email to the top of collection metadata list. */
  if (contactName && contactEmail) {
    collectionMetadataLinks.push({
      label: "Contact",
      testId: "collection-contact",
      url: `mailto:${contactEmail}`,
      value: contactName,
    });
  }

  /* Add collection links to collection metadata list. */
  /* Collection metadata reordered for preferred order of display. */
  const orderedLinks = sortCollectionLinks(links);
  for (const orderedLink of orderedLinks) {
    /* Build and add any valid collection metadata. */
    const collectionMetadataLink = buildCollectionMetadataLink(orderedLink);
    if (!collectionMetadataLink) continue;
    collectionMetadataLinks.push(collectionMetadataLink);
  }

  /* Return valid collection metadata. */
  return collectionMetadataLinks;
}

/**
 * Build display model of DOI link associated with a collection, if any. Display publication metadata if it has been
 * retrieved for the DOI, otherwise display the DOI link as is.
 * @params links - Links associated with a collection.
 * @returns Display model of DOI link.
 */
function buildDoiMetadataLink(
  doiLink: Link,
  summaryCitation?: string
): CollectionMetadataLink | undefined {
  // Build display model of DOI link.
  const doiMetadataLink = buildCollectionMetadataLink(doiLink);
  if (!doiMetadataLink) {
    return;
  }

  // If there's no summary citation for the collection, return the DOI link as is..
  if (!summaryCitation) {
    return doiMetadataLink;
  }

  // There's a summary citation link for the collection, update DOI link display.
  return {
    ...doiMetadataLink,
    label: "Publication DOI",
    value: summaryCitation,
  };
}

/**
 * @deprecated - supersede by buildCollectionMetadata once filter feature flag is removed (#1718).
 * @param links
 */
export function renderLinks(links: Link[]): (JSX.Element | null)[] {
  /* Reorder links into preferred order of display. */
  const orderedLinks = sortCollectionLinks(links);

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

/**
 * @deprecated - supersede by buildCollectionMetadata once filter feature flag is removed (#1718).
 * @param contact_name
 * @param contact_email
 */
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
      const assets = dataset.dataset_assets;
      const numOfDeployments = dataset.dataset_deployments.length;
      // Assets must contain a cxg and an h5ad. RDS (Seurat) are no longer mandatory for publishing
      return (
        numOfDeployments === 1 &&
        assets.some((asset) => asset.filetype === DATASET_ASSET_FORMAT.CXG) &&
        assets.some((asset) => asset.filetype === DATASET_ASSET_FORMAT.H5AD)
      );
    }) &&
    // (ebezzi): We need to ensure all dataset `processing_status` to be success, since creating a revision relies
    // on this status. Otherwise in the case of RDS is still WIP and a user publishes the collection and immediately
    // tries to create a revision, the revision creation will fail
    datasets.every(
      (dataset) =>
        dataset.processing_status.processing_status ===
        PROCESSING_STATUS.SUCCESS
    ) &&
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

/**
 * Returns collection metadata link or undefined if the metadata is incomplete.
 * @param link - link associated with collection.
 * @returns CollectionMetadataLink or undefined.
 */
function buildCollectionMetadataLink(
  link: Link
): CollectionMetadataLink | undefined {
  const { link_url: url, link_type: type, link_name: name } = link;

  /* Early exit; collection metadata link is incomplete when there is no corresponding link type. */
  const linkTypeOption = COLLECTION_LINK_TYPE_OPTIONS[type];
  if (!linkTypeOption) return;

  /* Grab the metadata display value. */
  let value: string | null = name;

  if (!value) {
    if (linkTypeOption === COLLECTION_LINK_TYPE_OPTIONS["DOI"]) {
      value = new URL(url).pathname.substring(1);
    } else {
      value = getUrlHost(url);
    }
  }

  /* Early exit; collection metadata link is incomplete when there is no display value. */
  if (!value) return;

  return {
    label: linkTypeOption.text,
    testId: "collection-link",
    url: url,
    value: value,
  };
}

function isPrivateRevision(collection: Collection) {
  return (
    collection.visibility === VISIBILITY_TYPE.PRIVATE && collection.has_revision
  );
}

/**
 * Sort collection links in preferred order of display.
 * @param links - links associated with collection.
 * @returns Array of collection links in preferred order of display.
 */
function sortCollectionLinks(links: Link[]): Link[] {
  const linkTypesToLinks = createLinkTypesToLinks(links);

  const orderedLinks: Link[] = [];

  for (const linkType of LINK_ORDER) {
    const links = linkTypesToLinks[linkType];

    if (!links) continue;

    for (const link of links) {
      orderedLinks.push(link);
    }
  }

  return orderedLinks;
}
