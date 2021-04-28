import React, { FC, ReactChild } from "react";
import { ACCESS_TYPE, VISIBILITY_TYPE } from "src/common/entities";
import {
  CollectionResponse,
  RevisionResponse,
  REVISION_STATUS,
} from "src/common/queries/collections";
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

function generateRevisions(
  collections: CollectionResponse[],
  revisionsEnabled: boolean
): RevisionResponse[] {
  const newCollections = collections.map(
    (collection): RevisionResponse => {
      return { ...collection, revision: REVISION_STATUS.DISABLED };
    }
  );
  // If revisions are disabled just return the array with objects' revisions status set to disabled
  if (!revisionsEnabled) return newCollections;

  const revisionMap = new Map<string, RevisionResponse>();

  newCollections.forEach((collection) => {
    const revisionObj = revisionMap.get(collection.id) || collection;
    if (revisionObj.revision !== REVISION_STATUS.STARTED) {
      revisionObj.revision =
        revisionObj.revision === REVISION_STATUS.DISABLED
          ? REVISION_STATUS.NOT_STARTED
          : REVISION_STATUS.STARTED;
      revisionMap.set(collection.id, revisionObj);
    }
  });
  return Array.from(revisionMap.values());
}

export default CollectionsGrid;
