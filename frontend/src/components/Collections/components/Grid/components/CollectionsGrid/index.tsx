import { FC, ReactChild } from "react";
import { ACCESS_TYPE, VISIBILITY_TYPE } from "src/common/entities";
import {
  CollectionResponse,
  CollectionResponsesMap,
} from "src/common/queries/collections";
import {
  CollectionHeaderCell,
  LeftAlignedHeaderCell,
  RightAlignedHeaderCell,
  StyledCollectionsGrid,
} from "../../common/style";
import CollectionRow from "../Row/CollectionRow";

interface Props {
  collections: CollectionResponsesMap;
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
  collections: CollectionResponsesMap,
  displayVisibility?: VISIBILITY_TYPE,
  accessType?: ACCESS_TYPE,
  revisionsEnabled = false
) {
  const collectionElements = [] as Array<ReactChild>;
  sortByCreatedAtDescendingAndFormat(collections).forEach(
    ({ id, visibility }) => {
      if (!displayVisibility || visibility === displayVisibility) {
        collectionElements.push(
          <CollectionRow
            id={id}
            key={id + visibility}
            visibility={visibility}
            accessType={accessType}
            revisionsEnabled={revisionsEnabled}
          />
        );
      }
    }
  );

  return collectionElements;
}

function sortByCreatedAtDescendingAndFormat(
  collectionsMap: CollectionResponsesMap
): CollectionResponse[] {
  return (
    Array.from(collectionsMap.values())
      .map((collectionsWithIDMap) => {
        return (collectionsWithIDMap.get(VISIBILITY_TYPE.PUBLIC) ||
          collectionsWithIDMap.get(
            VISIBILITY_TYPE.PRIVATE
          )) as CollectionResponse;
      })
      .sort((a, b) => b.created_at - a.created_at) || []
  );
}

export default CollectionsGrid;
