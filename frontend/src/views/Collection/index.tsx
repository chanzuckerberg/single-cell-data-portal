import { H3, Tab, Tabs } from "@blueprintjs/core";
import { RouteComponentProps } from "@reach/router";
import React, { FC, useState } from "react";
import { ACCESS_TYPE, VISIBILITY_TYPE } from "src/common/entities";
import { useCollection } from "src/common/queries/collections";
import DeleteCollection from "src/components/Collections/components/DeleteCollection";
import PublishCollection from "src/components/Collections/components/PublishCollection";
import { ViewGrid } from "../globalStyle";
import DatasetTab from "./components/DatasetTab";
import {
  CollectionButtons,
  CollectionInfo,
  Description,
  LinkContainer,
  TabWrapper,
} from "./style";
import { getIsPublishable, renderLinks } from "./utils";

interface RouteProps {
  id?: string;
}

export type Props = RouteComponentProps<RouteProps>;

const Collection: FC<Props> = ({ id = "" }) => {
  const isPrivate = window.location.pathname.includes("/private");

  const visibility = isPrivate
    ? VISIBILITY_TYPE.PRIVATE
    : VISIBILITY_TYPE.PUBLIC;

  const { data: collection, isError } = useCollection(id, visibility);

  const [selectedTab, setSelectedTab] = useState("datasets-tab");

  if (!collection || isError) {
    return null;
  }

  const datasets = collection.datasets;

  const isPublishable = getIsPublishable(datasets);

  const handleOnChange = function (newTabId: string) {
    setSelectedTab(newTabId);
  };

  return (
    <ViewGrid>
      <CollectionInfo>
        <H3>{collection.name}</H3>
        <Description>{collection.description}</Description>
        <LinkContainer>{renderLinks(collection.links)}</LinkContainer>
      </CollectionInfo>

      <CollectionButtons>
        {collection.access_type === ACCESS_TYPE.WRITE && isPrivate && (
          <DeleteCollection id={id} />
        )}
        {isPrivate && (
          <PublishCollection isPublishable={isPublishable} id={id} />
        )}
      </CollectionButtons>
      <TabWrapper>
        <Tabs
          onChange={handleOnChange}
          selectedTabId={selectedTab}
          id="collection-tabs"
          defaultSelectedTabId={"datasets-tab"}
        >
          <Tab
            id="datasets-tab"
            title="Datasets"
            panel={
              <DatasetTab
                collectionID={id}
                datasets={datasets}
                visibility={visibility}
              />
            }
          />
          <Tab id="geneset-tab" title="Genesets" panel={<div>HELLO</div>} />
        </Tabs>
      </TabWrapper>
    </ViewGrid>
  );
};

export default Collection;
