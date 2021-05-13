import React, { FC, ReactChild } from "react";
import { ACCESS_TYPE, VISIBILITY_TYPE } from "src/common/entities";
import {
  CollectionResponse,
  RevisionResponse,
} from "src/common/queries/collections";
import { generateRevisions } from "src/components/Collections/util";
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
  revisionsEnabled?: boolean;
}

const CollectionsGrid: FC<Props> = ({
  collections,
  accessType,
  displayVisibility,
  revisionsEnabled,
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
          {revisionsEnabled && (
            <RightAlignedHeaderCell></RightAlignedHeaderCell>
          )}
        </tr>
      </thead>
      <tbody>
        {renderCollections(
          collections,
          displayVisibility,
          accessType,
          revisionsEnabled
        )}
      </tbody>
    </StyledCollectionsGrid>
  );
};

function renderCollections(
  collections: CollectionResponse[],
  displayVisibility?: VISIBILITY_TYPE,
  accessType?: ACCESS_TYPE,
  revisionsEnabled = false
) {
  const collectionElements = [] as Array<ReactChild>;
  const revisions = generateRevisions(collections, revisionsEnabled);
  sortByCreatedAtDescending(revisions).forEach(
    ({ id, visibility, revision }) => {
      if (!displayVisibility || visibility === displayVisibility) {
        collectionElements.push(
          <CollectionRow
            id={id}
            key={id + visibility}
            visibility={visibility}
            accessType={accessType}
            revisionStatus={revision}
          />
        );
      }
    }
  );

  return collectionElements;
}

function sortByCreatedAtDescending(
  collections: RevisionResponse[]
): RevisionResponse[] {
  return collections?.sort((a, b) => b.created_at - a.created_at) || [];
}

export default CollectionsGrid;
