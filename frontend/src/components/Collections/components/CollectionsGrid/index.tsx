import React, { FC } from "react";
import { CollectionResponse } from "src/common/queries/collections";
import CollectionRow from "./CollectionRow";
import {
  CollectionHeaderCell,
  LeftAlignedHeaderCell,
  RightAlignedHeaderCell,
  StyledCollectionsGrid,
} from "./style";

interface Props {
  collections: CollectionResponse[];
  showStatus: boolean;
  showAllWritable: boolean;
}

const CollectionsGrid: FC<Props> = ({
  collections,
  showStatus,
  showAllWritable,
}) => {
  return (
    <StyledCollectionsGrid bordered>
      <thead>
        <tr>
          <CollectionHeaderCell>Collection</CollectionHeaderCell>
          <LeftAlignedHeaderCell>Organ</LeftAlignedHeaderCell>
          <LeftAlignedHeaderCell>Assay</LeftAlignedHeaderCell>
          <LeftAlignedHeaderCell>Species</LeftAlignedHeaderCell>
          <RightAlignedHeaderCell>Cell Count</RightAlignedHeaderCell>
          {showStatus && (
            <RightAlignedHeaderCell>Status</RightAlignedHeaderCell>
          )}
        </tr>
      </thead>
      <tbody>
        {collections?.map((collection) => (
          <CollectionRow
            id={collection.id}
            key={collection.id}
            showStatus={showStatus}
            showAllWritable={showAllWritable}
          />
        ))}
      </tbody>
    </StyledCollectionsGrid>
  );
};

export default CollectionsGrid;
