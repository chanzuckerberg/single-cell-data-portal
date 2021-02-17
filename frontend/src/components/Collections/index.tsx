import React, { FC } from "react";
import { VISIBILITY_TYPE } from "src/common/entities";
import { useCollections } from "src/common/queries/collections";
import CreateCollection from "../CreateCollectionModal";
import CollectionsGrid from "./components/Grid/components/CollectionsGrid";
import { TitleAndDescription, TitleWrapper } from "./style";

const Collections: FC = () => {
  const { isFetching, data: collections } = useCollections();

  if (isFetching && !collections) return <div>Loading collections...</div>;

  if (!collections) return <div>Sorry, we could not find any collections</div>;

  return (
    <>
      <TitleWrapper>
        <TitleAndDescription>
          <h1>Collections</h1>
          <p>Explore public collections of datasets or create your own.</p>
        </TitleAndDescription>
        <CreateCollection />
      </TitleWrapper>
      <CollectionsGrid
        collections={collections}
        displayVisibility={VISIBILITY_TYPE.PUBLIC}
      />
    </>
  );
};

export default Collections;
