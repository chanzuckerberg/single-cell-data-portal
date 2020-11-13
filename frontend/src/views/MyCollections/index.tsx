import { H1, Text } from "@blueprintjs/core";
import React, { FC } from "react";
import { useCollections } from "src/common/queries/collections";
import CollectionRow from "src/components/collections/components/CollectionRow";
import { ViewGrid } from "../globalStyle";
import {
  CollectionHeaderCell,
  CollectionsGrid,
  LeftAlignedHeaderCell,
  RightAlignedHeaderCell,
  StyledCreateCollection,
  TitleAndDescription,
} from "./style";

const MyCollections: FC = () => {
  const { data: collections, isFetching } = useCollections();
  if (isFetching && !collections) return null;
  if (!collections) return null;

  return (
    <ViewGrid>
      <TitleAndDescription>
        <H1>My Collections</H1>
        <Text>
          A list of collections you have either created, or have been shared
          with you
        </Text>
      </TitleAndDescription>
      <StyledCreateCollection />
      <CollectionsGrid bordered>
        <thead>
          <tr>
            <CollectionHeaderCell>Collection</CollectionHeaderCell>
            <LeftAlignedHeaderCell>Organ</LeftAlignedHeaderCell>
            <LeftAlignedHeaderCell>Assay</LeftAlignedHeaderCell>
            <LeftAlignedHeaderCell>Species</LeftAlignedHeaderCell>
            <RightAlignedHeaderCell>Cell Count</RightAlignedHeaderCell>
            <RightAlignedHeaderCell>Status</RightAlignedHeaderCell>
          </tr>
        </thead>
        <tbody>
          {collections?.map((collection) => (
            <CollectionRow id={collection.id} key={collection.id} />
          ))}
        </tbody>
      </CollectionsGrid>
    </ViewGrid>
  );
};

export default MyCollections;
