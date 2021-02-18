import React, { FC, ReactChild } from "react";
import { ACCESS_TYPE, VISIBILITY_TYPE } from "src/common/entities";
import { CollectionResponse } from "src/common/queries/collections";
import {
  CollectionHeaderCell,
  LeftAlignedHeaderCell,
  RightAlignedHeaderCell,
  StyledCollectionsGrid,
} from "../../common/style";
import CollectionRow from "../Row/CollectionRow";

interface Props {
  collections: CollectionResponse[];
  accessType?: ACCESS_TYPE;
  displayVisibility?: VISIBILITY_TYPE;
}

const CollectionsGrid: FC<Props> = ({
  collections,
  accessType,
  displayVisibility,
}) => {
  return (
    <StyledCollectionsGrid bordered>
      <thead>
        <tr>
          <CollectionHeaderCell>Collection</CollectionHeaderCell>
          <LeftAlignedHeaderCell>Tissue</LeftAlignedHeaderCell>
          <LeftAlignedHeaderCell>Assay</LeftAlignedHeaderCell>
          <LeftAlignedHeaderCell>Disease</LeftAlignedHeaderCell>
          <LeftAlignedHeaderCell>Organism</LeftAlignedHeaderCell>
          <RightAlignedHeaderCell>Cell Count</RightAlignedHeaderCell>
        </tr>
      </thead>
      <tbody>
        {sortByUpdatedAtDescending(collections)?.reduce(
          (acc, { id, visibility }) => {
            if (!displayVisibility || visibility === displayVisibility) {
              acc.push(
                <CollectionRow
                  id={id}
                  key={id + visibility}
                  visibility={visibility}
                  accessType={accessType}
                />
              );
            }
            return acc;
          },
          [] as Array<ReactChild>
        )}
      </tbody>
    </StyledCollectionsGrid>
  );
};

function sortByUpdatedAtDescending(
  collections: CollectionResponse[]
): CollectionResponse[] {
  return collections?.sort((a, b) => b.created_at - a.created_at) || [];
}

export default CollectionsGrid;
