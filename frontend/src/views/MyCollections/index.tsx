import { H1, HTMLTable, Text } from "@blueprintjs/core";
import React, { FC } from "react";
import { useCollections } from "src/common/queries/collections";
import CollectionRow from "src/components/collections/components/CollectionRow";
import { ViewGrid } from "../globalStyle";
import { StyledCreateCollection, TitleAndDescription } from "./style";

const MyCollections: FC<Props> = () => {
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
      <HTMLTable bordered style={{ gridColumn: "1/9" }}>
        <thead>
          <tr>
            <th>Collection</th>
            <th>Organ</th>
            <th>Assay</th>
            <th>Species</th>
            <th>Cell Count</th>
            <th>Status</th>
          </tr>
        </thead>
        <tbody>
          {collections?.map((collection) => (
            <CollectionRow id={collection.id} key={collection.id} />
          ))}
        </tbody>
      </HTMLTable>
    </ViewGrid>
  );
};

export default MyCollections;
