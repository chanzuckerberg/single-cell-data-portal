import { H3, Tab, Tabs } from "@blueprintjs/core";
import { RouteComponentProps } from "@reach/router";
import React, { FC, useState } from "react";
import { ACCESS_TYPE, VISIBILITY_TYPE } from "src/common/entities";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import { useCollection } from "src/common/queries/collections";
import DeleteCollection from "src/components/Collections/components/DeleteCollection";
import PublishCollection from "src/components/Collections/components/PublishCollection";
import DatasetTab from "src/views/Collection/components/DatasetTab";
import { ViewGrid } from "../globalStyle";
import GeneSetTab from "./components/GeneSetTab";
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

enum TABS {
  GENE_SETS = "gene-sets-tab",
  DATASETS = "datasets-tab",
}

export type Props = RouteComponentProps<RouteProps>;

const Collection: FC<Props> = ({ id = "" }) => {
  const isPrivate = window.location.pathname.includes("/private");

  const visibility = isPrivate
    ? VISIBILITY_TYPE.PRIVATE
    : VISIBILITY_TYPE.PUBLIC;

  const { data: collection, isError } = useCollection(id, visibility);

  const [selectedTab, setSelectedTab] = useState(TABS.DATASETS);

  if (!collection || isError) {
    return null;
  }

  const datasets = collection.datasets;

  const isPublishable = getIsPublishable(datasets);

  const handleOnChange = function (newTabId: TABS) {
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
          defaultSelectedTabId={TABS.DATASETS}
        >
          <Tab
            id={TABS.DATASETS}
            title="Datasets"
            panel={
              <DatasetTab
                collectionID={id}
                datasets={datasets}
                visibility={visibility}
              />
            }
          />
          {get(FEATURES.GENE_SETS) === BOOLEAN.TRUE && (
            <Tab id={TABS.GENE_SETS} title="Gene Sets" panel={<GeneSetTab />} />
          )}
        </Tabs>
      </TabWrapper>
    </ViewGrid>
  );
};

export default Collection;
