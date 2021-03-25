import { Intent, Tag } from "@blueprintjs/core";
import loadable from "@loadable/component";
import Link from "next/link";
import React, { FC } from "react";
import {
  ACCESS_TYPE,
  COLLECTION_LINK_TYPE,
  VISIBILITY_TYPE,
} from "src/common/entities";
import { useCollection } from "src/common/queries/collections";
import { aggregateDatasetsMetadata } from "../../../common/utils";
import {
  LeftAlignedDetailsCell,
  RightAlignedDetailsCell,
  StyledCell,
  StyledRow,
} from "../common/style";
import { CollectionTitleText, ContactText, DOILink } from "./style";

interface Props {
  id: string;
  accessType?: ACCESS_TYPE;
  visibility: VISIBILITY_TYPE;
}

const AsyncPopover = loadable(
  () =>
    /*webpackChunkName: 'CollectionRow/components/Popover' */ import(
      "src/components/Collections/components/Grid/components/Popover"
    )
);

const conditionalPopover = (values: string[]) => {
  if (!values || values.length === 0) {
    return <LeftAlignedDetailsCell>-</LeftAlignedDetailsCell>;
  }

  return <AsyncPopover values={values} />;
};

type DOI = {
  doi: string;
  link: string;
};

const CollectionRow: FC<Props> = (props) => {
  const { data: collection } = useCollection({
    id: props.id,
    visibility: props.visibility,
  });

  if (!collection) return null;

  // If there is an explicity accessType only show collections with that accessType
  if (props.accessType && collection.access_type !== props.accessType) {
    return null;
  }

  const { id, links, visibility, contact_name, name, datasets } = collection;

  const dois: Array<DOI> = links.reduce((acc, link) => {
    if (link.link_type !== COLLECTION_LINK_TYPE.DOI) return acc;

    const url = new URL(link.link_url);

    acc.push({ doi: url.pathname.substring(1), link: link.link_url });

    return acc;
  }, [] as DOI[]);

  const isPrivate = visibility === VISIBILITY_TYPE.PRIVATE;

  const {
    tissue,
    assay,
    disease,
    organism,
    cell_count,
  } = aggregateDatasetsMetadata(datasets);

  return (
    <StyledRow>
      <StyledCell>
        <Link href={`/collections/${id}${isPrivate ? "/private" : ""}`}>
          <CollectionTitleText
            href={`/collections/${id}${isPrivate ? "/private" : ""}`}
            data-test-id="collection-link"
          >
            {name}
          </CollectionTitleText>
        </Link>
        <ContactText>{contact_name}</ContactText>

        {props.accessType === ACCESS_TYPE.WRITE ? (
          <Tag
            minimal
            intent={isPrivate ? Intent.PRIMARY : Intent.SUCCESS}
            data-test-id="visibility-tag"
          >
            {isPrivate ? "Private" : "Published"}
          </Tag>
        ) : (
          dois?.map((doi) => (
            <DOILink
              key={doi.doi}
              href={doi.link}
              target="_blank"
              rel="noopener"
            >
              {doi.doi}
            </DOILink>
          ))
        )}
      </StyledCell>
      {conditionalPopover(tissue)}
      {conditionalPopover(assay)}
      {conditionalPopover(disease)}
      {conditionalPopover(organism)}
      <RightAlignedDetailsCell>{cell_count || "-"}</RightAlignedDetailsCell>
    </StyledRow>
  );
};

export default CollectionRow;
