import { H1, Text } from "@blueprintjs/core";
import React, { FC } from "react";
import { useCollections } from "src/common/queries/collections";
import { ViewGrid } from "../globalStyle";
import { StyledCreateCollection, TitleAndDescription } from "./style";

const MyCollections: FC<Props> = () => {
  const { data: collection, isLoading } = useCollections();
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
    </ViewGrid>
  );
};

export default MyCollections;
