import React, { FC } from "react";
import { VISIBILITY_TYPE } from "src/common/entities";
import { useCollections } from "src/common/queries/collections";
import CreateCollection from "../CreateCollectionModal";
import CollectionsGrid from "./components/Grid/components/CollectionsGrid";
import { TitleWrapper } from "./style";

const Collections: FC = () => {
  const { isFetching, data: collections } = useCollections();

  if (isFetching && !collections) return <div>Loading collections...</div>;

  if (!collections) return <div>Sorry, we could not find any collections</div>;

  return (
    <>
      <TitleWrapper>
        <h1>Collections</h1>
        <CreateCollection />
      </TitleWrapper>
      <p>Explore public collections of datasets or create your own.</p>
      <CollectionsGrid
        collections={collections}
        displayVisibility={VISIBILITY_TYPE.PUBLIC}
      />
    </>
  );
};

export default Collections;
