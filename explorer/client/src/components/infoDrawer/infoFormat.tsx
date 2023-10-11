// Core dependencies
import { H3, HTMLTable, Classes } from "@blueprintjs/core";
import React, { CSSProperties } from "react";

// App dependencies
import {
  Author,
  Consortium,
  DatasetMetadata,
  Link,
  PublisherMetadata,
} from "../../common/types/entities";
import { Category } from "../../common/types/schema";
import * as globals from "../../globals";

const COLLECTION_LINK_ORDER_BY = [
  "DOI",
  "DATA_SOURCE",
  "RAW_DATA",
  "PROTOCOL",
  "LAB_WEBSITE",
  "OTHER",
];

interface LinkView {
  name: string;
  type: string;
  url: string;
}

interface MetadataView {
  key: string;
  value: string;
}

interface Props {
  datasetMetadata: DatasetMetadata;
  singleValueCategories: SingleValueCategories;
}

export type SingleValueCategories = Map<string, Category>;

/**
 * Sort collection links by custom sort order, create view-friendly model of link types.
 * @param links - Links associated with a collection.
 * @param contactName - Contact name associated with a collection.
 * @param contactEmail - Contact email associated with a collection.
 * @param publisherMetadata - Publication metadata associated with a collection.
 * @returns Array of link objects formatted for display.
 */
const buildCollectionLinks = (
  links: Link[],
  contactName: string,
  contactEmail: string,
  publisherMetadata?: PublisherMetadata,
): LinkView[] => {
  const collectionMetadataLinks = [];

  // Clone the links so we can safely remove the DOI link type if present and display it separately from the other link
  // types at the top of the metadata links list.
  const linksClone = [...links];

  // If collection has an associated DOI, display either the summary citation or the DOI itself.
  const doiLink = links.find((link: Link) => link.link_type === "DOI");
  if (doiLink) {
    const doiMetadataLink = buildDoiMetadataLink(doiLink, publisherMetadata);
    collectionMetadataLinks.push(doiMetadataLink);

    // Remove the DOI from links so we don't display it again in the links section below contact.
    linksClone.splice(linksClone.indexOf(doiLink), 1);
  }

  // Add contact name and email to the top of collection metadata list.
  if (contactName && contactEmail) {
    collectionMetadataLinks.push({
      type: "Contact",
      url: `mailto:${contactEmail}`,
      name: contactName,
    });
  }

  const sortedLinks = [...linksClone].sort(sortCollectionLinks);
  sortedLinks.forEach((link) => {
    const { link_name: name, link_type: type, link_url: url } = link;
    collectionMetadataLinks.push({
      name: buildLinkName(name, type, url),
      type: transformLinkTypeToDisplay(type),
      url,
    });
  });

  return collectionMetadataLinks;
};

/**
 * Build display model of DOI link associated with a collection, if any. Display publication metadata if it has been
 * retrieved for the DOI, otherwise display the DOI link as is.
 * @param doiLink - Link with type DOI, associated with a collection.
 * @param publisherMetadata - Publication metadata associated with collection.
 * @returns Display model of DOI link.
 */
const buildDoiMetadataLink = (
  doiLink: Link,
  publisherMetadata?: PublisherMetadata,
): LinkView => {
  // Build display model of DOI link.
  const { link_name: name, link_type: type, link_url: url } = doiLink;
  const doiMetadataLink = {
    name: buildLinkName(name, type, url),
    type: transformLinkTypeToDisplay(type),
    url,
  };

  // If there's no summary citation for the collection, return the DOI link as is.
  const summaryCitation = buildSummaryCitation(publisherMetadata);
  if (!summaryCitation) {
    return doiMetadataLink;
  }

  // There's a summary citation link for the collection; update DOI link display.
  return {
    ...doiMetadataLink,
    type: "Publication",
    name: summaryCitation,
  };
};

/**
 * Determine name to display for collection link.
 * @param name - Link display name
 * @param type - Link type (e.g. DOI)
 * @param url - Link URL
 * @returns Pathname if link type is DOI otherwise host.
 */
const buildLinkName = (name: string, type: string, url: string): string => {
  if (name) {
    return name;
  }
  let validUrl;
  try {
    validUrl = new URL(url);
  } catch (e) {
    return url;
  }
  if (type === "DOI") {
    return validUrl.pathname.substring(1);
  }
  return validUrl.host;
};

/**
 * Build summary citation format from given publisher metadata:
 * Last name of first author (publication year) journal abbreviation such as Ren et al. (2021) Cell.
 * @param publisherMetadata - Publication metadata of a collection.
 */
const buildSummaryCitation = (
  publisherMetadata?: PublisherMetadata,
): string => {
  if (!publisherMetadata) {
    return "";
  }

  const citationTokens = [];

  // Add author to citation - family name if first author is a person, name if first author is a consortium.
  const { authors, journal, published_year: publishedYear } = publisherMetadata;
  const [firstAuthor] = authors;
  if (firstAuthor) {
    if (isAuthorPerson(firstAuthor)) {
      citationTokens.push(firstAuthor.family);
    } else {
      citationTokens.push(firstAuthor.name);
    }

    if (authors.length > 1) {
      citationTokens.push("et al.");
    }
  }

  // Add year and journal.
  citationTokens.push(`(${publishedYear})`);
  citationTokens.push(journal);

  return citationTokens.join(" ");
};

/**
 * Generate inline styles to be applied to collections and meta tables.
 * @returns Inline style object.
 */
const getTableStyles = (): CSSProperties => ({
  tableLayout: "fixed",
  width: "100%",
});

/**
 * Publication authors can be a person or a consortium; determine if the given author is in fact a person (and not a
 * consortium).
 * @param author - Person or consortium associated with a publication.
 * @returns True if author is a person and not a consortium.
 */
const isAuthorPerson = (author: Author | Consortium): author is Author =>
  (author as Author).given !== undefined;

/**
 * Render collection contact and links.
 * @param datasetMetadata - Dataset metadata containing collection link information to be displayed
 * @returns Markup displaying contact and collection-related links.
 */
const renderCollectionLinks = (
  datasetMetadata: DatasetMetadata,
): JSX.Element => {
  const {
    collection_contact_name: contactName,
    collection_contact_email: contactEmail,
  } = datasetMetadata;
  const links = buildCollectionLinks(
    datasetMetadata.collection_links,
    contactName,
    contactEmail,
    datasetMetadata.collection_publisher_metadata,
  );

  return (
    <>
      {renderSectionTitle("Collection")}
      <HTMLTable style={getTableStyles()}>
        <tbody>
          {links.map(({ name, type, url }, i) => (
            <tr {...{ key: i }}>
              <td>{type}</td>
              <td>
                <a href={url} rel="noopener" target="_blank">
                  {name}
                </a>
              </td>
            </tr>
          ))}
        </tbody>
      </HTMLTable>
    </>
  );
};

/**
 * Render dataset metadata. That is, attributes found in categorical fields.
 * @param singleValueCategories - Attributes from categorical fields
 * @returns Markup for displaying meta in table format.
 */
const renderDatasetMetadata = (
  singleValueCategories: SingleValueCategories,
): JSX.Element | null => {
  if (singleValueCategories.size === 0) {
    return null;
  }
  const metadataViews = buildDatasetMetadataViews(singleValueCategories);
  metadataViews.sort(sortDatasetMetadata);
  return (
    <>
      {renderSectionTitle("Dataset")}
      <HTMLTable style={getTableStyles()}>
        <tbody>
          {metadataViews.map(({ key, value }) => (
            <tr {...{ key }}>
              <td>{key}</td>
              <td>{value}</td>
            </tr>
          ))}
        </tbody>
      </HTMLTable>
    </>
  );
};

/**
 * Create DOM elements for displaying section title.
 * @param title - Section title to display
 * @returns Styled markup representation displaying section title.
 */
const renderSectionTitle = (title: string): JSX.Element => (
  <p style={{ margin: "24px 0 8px" }}>
    <strong>{title}</strong>
  </p>
);

/**
 * Compare function for sorting collection links by custom link type order.
 * @param link0 - First link to compare
 * @param link1 - Second link value to compare
 * @returns Number indicating sort precedence of link0 vs link1.
 */
const sortCollectionLinks = (link0: Link, link1: Link): number =>
  COLLECTION_LINK_ORDER_BY.indexOf(link0.link_type) -
  COLLECTION_LINK_ORDER_BY.indexOf(link1.link_type);

/**
 * Compare function for metadata key value pairs by key - alpha, ascending.
 * @param m0 - First metadata value to compare
 * @param m1 - Second metadata value to compare
 * @returns Number indicating sort precedence of m0 vs m1.
 */
const sortDatasetMetadata = (m0: MetadataView, m1: MetadataView) => {
  if (m0.key < m1.key) {
    return -1;
  }
  if (m0.key > m1.key) {
    return 1;
  }
  return 0;
};

/**
 * Convert link type from upper snake case to title case.
 * @param type - Link type to transform.
 * @returns Transformed link type ready for display.
 */
const transformLinkTypeToDisplay = (type: string): string => {
  const tokens = type.split("_");
  return tokens
    .map((token) => token.charAt(0) + token.slice(1).toLowerCase())
    .join(" ");
};

/**
 * Build array of view model objects from given single value categories map, ignoring ontology terms or metadata
 * without values. Add ontology terms as tooltips of their corresponding values.
 * @param singleValueCategories - Attributes from categorical fields
 * @returns  Array of metadata key/value pairs.
 */
const buildDatasetMetadataViews = (
  singleValueCategories: SingleValueCategories,
): MetadataView[] =>
  Array.from(singleValueCategories.entries())
    .filter(([key, value]) => {
      if (key.indexOf(globals.ONTOLOGY_KEY) >= 0) {
        // skip ontology terms
        return false;
      }
      // skip metadata without values
      return value;
    })
    .map(([key, value]) => ({ key, value: String(value) }));

const InfoFormat = React.memo<Props>(
  ({ datasetMetadata, singleValueCategories }) => (
    <div className={Classes.DRAWER_BODY}>
      <div className={Classes.DIALOG_BODY}>
        <H3>{datasetMetadata.collection_name}</H3>
        <p>{datasetMetadata.collection_description}</p>
        {renderCollectionLinks(datasetMetadata)}
        {renderDatasetMetadata(singleValueCategories)}
      </div>
    </div>
  ),
);

export default InfoFormat;
