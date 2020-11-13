import { Button, Classes, Popover, Position } from "@blueprintjs/core";
import React, { FC } from "react";
import { COLLECTION_LINK_TYPE } from "src/common/entities";
import { useCollection } from "src/common/queries/collections";
import {
  CollectionTitleText,
  DetailsCell,
  DOIText,
  StyledCell,
  StyledRow,
} from "./style";

interface Props {
  id: string;
}

const conditionalPopover = (values: string[]) => {
  if (!values || values.length === 0) return "-";
  return (
    <>
      <p>
        {values[0]}
        <br />
        {values[1]}
      </p>
      {values.length > 2 && (
        <Popover
          boundary="window"
          modifiers={{
            preventOverflow: { enabled: false },
          }}
          usePortal
          position={Position.BOTTOM}
          popoverClassName={Classes.POPOVER_CONTENT_SIZING}
          content={
            <div style={{ margin: "16px 16px" }}>
              <p>
                {values.map((val, idx) => (
                  <>
                    {val}
                    {idx !== values.length - 1 && <br />}
                  </>
                ))}
              </p>
            </div>
          }
        >
          <Button minimal>+{values.length - 2}</Button>
        </Popover>
      )}
    </>
  );
};

const CollectionRow: FC<Props> = ({ id }) => {
  const { data: collection } = useCollection(id);

  if (!collection) return null;

  const dois = collection.links.reduce((acc, link) => {
    if (link.type !== COLLECTION_LINK_TYPE.DOI) return acc;
    const url = new URL(link.url);
    acc.push(url.pathname.substring(1));
    return acc;
  }, [] as string[]);

  /* DEMO */
  if (dois.length === 0) {
    dois.push("10.1126/science.aax6234");
  }
  return (
    <StyledRow>
      <StyledCell>
        <CollectionTitleText href={`/collections/${id}`}>
          {collection.name}
        </CollectionTitleText>
        <div>{collection.contact_name ?? "Seve Badajoz"}</div>
        {dois?.map((doi) => (
          <DOIText key={doi}>{doi}</DOIText>
        ))}
      </StyledCell>
      <DetailsCell>{conditionalPopover(["tummy", "bum"])}</DetailsCell>
      <DetailsCell>
        {conditionalPopover(["tummy", "bum", "bum", "bum"])}
      </DetailsCell>
      <DetailsCell>{conditionalPopover(["tummy"])}</DetailsCell>
      {/* <DetailsCell>{conditionalPopover(collection.assays)}</DetailsCell> */}

      <DetailsCell>{conditionalPopover(collection.species)}</DetailsCell>
      <DetailsCell>{collection.cell_count || "-"}</DetailsCell>
      <DetailsCell>{collection.status || "-"}</DetailsCell>
    </StyledRow>
  );
};

export default CollectionRow;
