import { H3, Intent, Tab, Tabs } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import Head from "next/head";
import { useRouter } from "next/router";
import { FC, useEffect, useState } from "react";
import { ROUTES } from "src/common/constants/routes";
import { ACCESS_TYPE, VISIBILITY_TYPE } from "src/common/entities";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { useExplainNewTab } from "src/common/hooks/useExplainNewTab";
import { BOOLEAN } from "src/common/localStorage/set";
import {
  useCollection,
  useCollectionUploadLinks,
} from "src/common/queries/collections";
import { UploadingFile } from "src/components/DropboxChooser";
import DatasetTab from "src/views/Collection/components/DatasetTab";
import { ViewGrid } from "../globalStyle";
import ActionButtons from "./components/ActionButtons";
import Banner from "./components/Banner";
import GeneSetTab from "./components/GeneSetTab";
import Toast from "./components/Toast";
import {
  CollectionInfo,
  Description,
  LinkContainer,
  StyledCallout,
  TabWrapper,
} from "./style";
import {
  getIsPublishable,
  renderContact,
  renderLinks,
  revisionIsPublishable,
} from "./utils";

enum TABS {
  GENE_SETS = "gene-sets-tab",
  DATASETS = "datasets-tab",
}

const Collection: FC = () => {
  const router = useRouter();
  const { params, tombstoned_dataset_id } = router.query;

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

  if (collection?.tombstone === true) {
    router.push(ROUTES.HOMEPAGE + "?tombstoned_dataset_id=" + id);
  }

  const isCurator = get(FEATURES.CURATOR) === BOOLEAN.TRUE;
  const isRevision = isCurator && !!collection?.has_revision;

  const [selectedTab, setSelectedTab] = useState(TABS.DATASETS);

  const collectionContactName = collection?.contact_name;
  useEffect(() => {
    // collectionContactName would be undefined in the case where collection is undefined.
    if (!tombstoned_dataset_id || typeof collectionContactName === "undefined")
      return;

    let message = "";
    // TODO: Remove empty string check after re-curation is complete. Contact name should always be populated.
    if (collectionContactName) {
      message = `A dataset was withdrawn by ${collectionContactName}. You've been redirected to the parent collection.`;
    } else {
      message =
        "A dataset was withdrawn. You've been redirected to the parent collection.";
    }

    Toast.show({
      icon: IconNames.ISSUE,
      intent: Intent.PRIMARY,
      message,
    });
  }, [tombstoned_dataset_id, collectionContactName]);

  /* Pop toast if user has come from Explorer with work in progress */
  useExplainNewTab(
    "To maintain your in-progress work on the previous dataset, we opened this collection in a new tab."
  );

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
          Toast.show({
            icon: IconNames.TICK,
            intent: Intent.PRIMARY,
            message:
              "Your file is being uploaded which will continue in the background, even if you close this window.",
          });
        },
      }
    );
  };

  let datasets = Array.from(collection.datasets.values());

  // Filter out tombstoned datasets if we're not looking at revision
  if (!isRevision) datasets = datasets.filter((dataset) => !dataset.tombstone);

  // (thuang): `isFetching` is used to prevent incorrectly enabling "Publish"
  // when React Query is fetching cached `collection` and its outdated
  // `datasets`
  const isPublishable =
    getIsPublishable(datasets) &&
    !isUploadingLink &&
    !isFetching &&
    revisionIsPublishable(collection, isCurator);

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
        {collection.has_revision && visibility === VISIBILITY_TYPE.PRIVATE && (
          <StyledCallout intent={Intent.PRIMARY} icon={null}>
            <span data-test-id="revision-status">
              {collection.revision_diff
                ? "This collection has changed since you last published it."
                : "This is a private revision of a public collection."}
            </span>
          </StyledCallout>
        )}

        {!isRevision && <Banner collectionId={collection.id} />}

        <CollectionInfo>
          <H3 data-test-id="collection-name">{collection.name}</H3>
          <Description data-test-id="collection-description">
            {collection.description}
          </Description>
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
            isRevision={isRevision}
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
                  isRevision={isRevision}
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
