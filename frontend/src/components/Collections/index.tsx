import { FC } from "react";
import { VISIBILITY_TYPE } from "src/common/entities";
import { useExplainNewTab } from "src/common/hooks/useExplainNewTab";
import { useCollections } from "src/common/queries/collections";
import CreateCollection from "../CreateCollectionModal";
import CollectionsGrid from "./components/Grid/components/CollectionsGrid";
import { TitleAndDescription, TitleWrapper } from "./style";

const Collections: FC = () => {
  const { isFetching, data: collections } = useCollections();

  /* Pop toast if user has come from Explorer with work in progress */
  useExplainNewTab(
    "To maintain your in-progress work on the previous dataset, we opened a new tab."
  );

  if (isFetching && !collections) return <div>Loading collections...</div>;

  if (!collections) return <div>Sorry, we could not find any collections</div>;

  return (
    <>
      <TitleWrapper>
        <TitleAndDescription>
          <h1 data-test-id="collections-header">Collections</h1>
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
