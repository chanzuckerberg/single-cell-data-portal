import { H1, Text } from "@blueprintjs/core";
import Head from "next/head";
import { FC } from "react";
import { ACCESS_TYPE } from "src/common/entities";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import { useCollections } from "src/common/queries/collections";
import CollectionsGrid from "src/components/Collections/components/Grid/components/CollectionsGrid";
import { ViewGrid } from "../globalStyle";
import { StyledCreateCollection, TitleAndDescription } from "./style";

const MyCollections: FC = () => {
  const { data: collections, isFetching } = useCollections();

  if (isFetching && !collections) return null;

  if (!collections) return null;

  // This prop is drilled down two levels to CollectionRow
  const revisionsEnabled = get(FEATURES.REVISION) === BOOLEAN.TRUE;

  return (
    <>
      <Head>
        <title>cellxgene | My Collections</title>
      </Head>
      <ViewGrid>
        <TitleAndDescription>
          <H1>My Collections</H1>
          <Text>A list of collections you have created</Text>
        </TitleAndDescription>
        <StyledCreateCollection />
        <CollectionsGrid
          collections={collections}
          accessType={ACCESS_TYPE.WRITE}
          revisionsEnabled={revisionsEnabled}
        />
      </ViewGrid>
    </>
  );
};

export default MyCollections;
