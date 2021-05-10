import { H3, Intent, Tab, Tabs } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import Head from "next/head";
import { useRouter } from "next/router";
import React, { FC, useState } from "react";
import { useQueryCache } from "react-query";
import { ACCESS_TYPE, VISIBILITY_TYPE } from "src/common/entities";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import {
  RevisionResponse,
  REVISION_STATUS,
  useCollection,
  useCollections,
  useCollectionUploadLinks,
  USE_COLLECTION,
} from "src/common/queries/collections";
import { generateRevisionMap } from "src/components/Collections/util";
import { UploadingFile } from "src/components/DropboxChooser";
import DatasetTab from "src/views/Collection/components/DatasetTab";
import { ViewGrid } from "../globalStyle";
import ActionButtons from "./components/ActionButtons";
import DatasetUploadToast from "./components/DatasetUploadToast";
import GeneSetTab from "./components/GeneSetTab";
import {
  CollectionInfo,
  Description,
  LinkContainer,
  TabWrapper,
} from "./style";
import { getIsPublishable, renderContact, renderLinks } from "./utils";

enum TABS {
  GENE_SETS = "gene-sets-tab",
  DATASETS = "datasets-tab",
}

const Collection: FC = () => {
  const router = useRouter();
  const { params } = router.query;

  let id = "";
  let isPrivate = false;

  if (Array.isArray(params)) {
    id = params[0];
    isPrivate = params[1] === "private";
  } else if (params) {
    id = params;
  }

  const visibility = isPrivate
    ? VISIBILITY_TYPE.PRIVATE
    : VISIBILITY_TYPE.PUBLIC;

  const [uploadLink] = useCollectionUploadLinks(id, visibility);

  const [isUploadingLink, setIsUploadingLink] = useState(false);

  const collectionState = useCollection({
    id,
    visibility,
  });

  const { data: collection, isError, isFetching } = collectionState;

  const { data: collections } = useCollections();
  const revisionsEnabled = get(FEATURES.REVISION) === BOOLEAN.TRUE;
  let isRevision = false;

  if (revisionsEnabled && isPrivate && collection && collections) {
    const revisionMap = generateRevisionMap(
      collections,
      revisionsEnabled
    ) as Map<string, RevisionResponse>;
    isRevision =
      revisionMap.get(collection.id)?.revision === REVISION_STATUS.STARTED;
  }

  const [selectedTab, setSelectedTab] = useState(TABS.DATASETS);

  const queryCache = useQueryCache();

  if (!collection || isError) {
    return null;
  }

  const addNewFile = (newFile: UploadingFile) => {
    if (!newFile.link) return;

    const payload = JSON.stringify({ url: newFile.link });

    setIsUploadingLink(true);

    uploadLink(
      { collectionId: id, payload },
      {
        onSettled: () => setIsUploadingLink(false),
        onSuccess: () => {
          DatasetUploadToast.show({
            icon: IconNames.TICK,
            intent: Intent.PRIMARY,
            message:
              "Your file is being uploaded which will continue in the background, even if you close this window.",
          });

          queryCache.invalidateQueries(USE_COLLECTION);
        },
      }
    );
  };

  const datasets = collection.datasets;

  // (thuang): `isFetching` is used to prevent incorrectly enabling "Publish"
  // when React Query is fetching cached `collection` and its outdated
  // `datasets`
  const isPublishable =
    getIsPublishable(datasets) && !isUploadingLink && !isFetching;

  const handleOnChange = function (newTabId: TABS) {
    setSelectedTab(newTabId);
  };

  const shouldShowPrivateWriteAction =
    collection.access_type === ACCESS_TYPE.WRITE && isPrivate;

  return (
    <>
      <Head>
        <title>cellxgene | {collection.name}</title>
      </Head>
      <ViewGrid>
        <CollectionInfo>
          <H3>{collection.name}</H3>
          <Description>{collection.description}</Description>
          <LinkContainer>
            {/*
             * (thuang): Contact and Links order defined here:
             * https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell/124
             */}
            {renderContact(collection.contact_name, collection.contact_email)}
            {renderLinks(collection.links)}
          </LinkContainer>
        </CollectionInfo>

        {shouldShowPrivateWriteAction && (
          <ActionButtons
            id={id}
            addNewFile={addNewFile}
            isPublishable={isPublishable}
          />
        )}

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
              <Tab
                id={TABS.GENE_SETS}
                title="Gene Sets"
                panel={<GeneSetTab />}
              />
            )}
          </Tabs>
        </TabWrapper>
      </ViewGrid>
    </>
  );
};

export default Collection;
