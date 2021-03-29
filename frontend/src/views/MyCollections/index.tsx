import { H1, Text } from "@blueprintjs/core";
import Head from "next/head";
import React, { FC } from "react";
import { ACCESS_TYPE } from "src/common/entities";
import { useCollections } from "src/common/queries/collections";
import CollectionsGrid from "src/components/Collections/components/Grid/components/CollectionsGrid";
import { ViewGrid } from "../globalStyle";
import { StyledCreateCollection, TitleAndDescription } from "./style";

const MyCollections: FC = () => {
  const { data: collections, isFetching } = useCollections();

  if (isFetching && !collections) return null;

  if (!collections) return null;

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
        />
      </ViewGrid>
    </>
  );
};

export default MyCollections;
