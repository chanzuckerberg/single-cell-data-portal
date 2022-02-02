import React from "react";
import {
  CollectionDescription as Description,
  CollectionDetails as Details,
} from "src/components/Collection/components/CollectionDetails/style";
import CollectionMetadata, {
  CollectionMetadataLink,
} from "src/components/Collection/components/CollectionMetadata";

interface Props {
  collectionDescription: string;
  collectionMetadata: CollectionMetadataLink[];
}

export default function CollectionDetails({
  collectionDescription,
  collectionMetadata,
}: Props): JSX.Element {
  return (
    <Details>
      {/* Collection Description*/}
      <Description data-test-id="collection-description">
        {collectionDescription}
      </Description>
      {/* Collection Metadata*/}
      {collectionMetadata && (
        <CollectionMetadata collectionMetadata={collectionMetadata} />
      )}
    </Details>
  );
}
