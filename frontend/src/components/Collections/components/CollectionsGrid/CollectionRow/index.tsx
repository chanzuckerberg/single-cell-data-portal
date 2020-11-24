import { Classes, Intent, Position, Tag } from "@blueprintjs/core";
import React, { FC } from "react";
import {
  ACCESS_TYPE,
  COLLECTION_LINK_TYPE,
  VISIBILITY_TYPE,
} from "src/common/entities";
import { useCollection } from "src/common/queries/collections";
import {
  CollectionTitleText,
  FieldValues,
  LeftAlignedDetailsCell,
  RightAlignedDetailsCell,
  StyledButton,
  StyledCell,
  StyledPopover,
  StyledRow,
} from "./style";

interface Props {
  id: string;
  accessType: ACCESS_TYPE;
  visibility: VISIBILITY_TYPE;
}

const conditionalPopover = (values: string[]) => {
  if (!values || values.length === 0) {
    return <LeftAlignedDetailsCell>-</LeftAlignedDetailsCell>;
  }

  return (
    <LeftAlignedDetailsCell>
      <FieldValues>
        {values[0]}
        <br />
        {values[1]}
      </FieldValues>
      {values.length > 2 && (
        <StyledPopover
          boundary="window"
          modifiers={{
            hide: { enabled: false },
            preventOverflow: { enabled: false },
          }}
          position={Position.BOTTOM}
          popoverClassName={Classes.POPOVER_CONTENT_SIZING}
          content={
            <div>
              {values.map((val, idx) => (
                <React.Fragment key={val}>
                  {val}
                  {idx !== values.length - 1 && <br />}
                </React.Fragment>
              ))}
            </div>
          }
        >
          <StyledButton minimal>+{values.length - 2}</StyledButton>
        </StyledPopover>
      )}
    </LeftAlignedDetailsCell>
  );
};

type DOI = {
  doi: string;
  link: string;
};

const CollectionRow: FC<Props> = ({ id, accessType, visibility }) => {
  const { data: collection } = useCollection(id, visibility);

  if (!collection) return null;

  // If there is an explicity accessType only show collections with that accessType
  if (accessType && collection.access_type !== accessType) return null;

  const dois: Array<DOI> = collection.links.reduce((acc, link) => {
    if (link.link_type !== COLLECTION_LINK_TYPE.DOI) return acc;

    const url = new URL(link.link_url);

    acc.push({ doi: url.pathname.substring(1), link: link.link_url });

    return acc;
  }, [] as DOI[]);

  const isPrivate = collection.visibility === VISIBILITY_TYPE.PRIVATE;

  // TODO: Generate data from datasets #737
  const { tissue, assays, disease, organism, cellCount } = {} as any;

  return (
    <StyledRow>
      <StyledCell>
        <CollectionTitleText
          href={`/collections/${id}${isPrivate ? "/private" : ""}`}
        >
          {collection.name}
        </CollectionTitleText>
        <div>{collection.contact_name}</div>
        {dois?.map((doi) => (
          <a key={doi.doi} href={doi.link}>
            {doi.doi}
          </a>
        ))}
        <Tag minimal intent={isPrivate ? Intent.PRIMARY : Intent.SUCCESS}>
          {isPrivate ? "Private" : "Published"}
        </Tag>
      </StyledCell>
      {conditionalPopover(tissue)}
      {conditionalPopover(assays)}
      {conditionalPopover(disease)}
      {conditionalPopover(organism)}
      <RightAlignedDetailsCell>{cellCount || "-"}</RightAlignedDetailsCell>
    </StyledRow>
  );
};

export default CollectionRow;
