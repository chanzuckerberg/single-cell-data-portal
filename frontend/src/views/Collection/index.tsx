import { Intent } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import Head from "next/head";
import { useRouter } from "next/router";
import React, { FC, useEffect, useState } from "react";
import { ROUTES } from "src/common/constants/routes";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { useExplainNewTab } from "src/common/hooks/useExplainNewTab";
import { BOOLEAN } from "src/common/localStorage/set";
import { useCollection } from "src/common/queries/collections";
import { removeParams } from "src/common/utils/removeParams";
import { isTombstonedCollection } from "src/common/utils/typeGuards";
import CollectionDescription from "src/components/Collection/components/CollectionDescription";
import CollectionMetadata from "src/components/Collection/components/CollectionMetadata";
import CollectionRevisionStatusCallout from "src/components/Collection/components/CollectionRevisionStatusCallout";
import DatasetTab from "src/views/Collection/components/DatasetTab";
import Toast from "./components/Toast";
import {
  CollectionConsortia,
  CollectionDetail,
  CollectionHero,
  CollectionView,
} from "./style";
import {
  buildCollectionMetadataLinks,
  getIsPublishable,
  isCollectionHasPrivateRevision,
  isCollectionPrivateRevision,
  revisionIsPublishable,
} from "./utils";
import CollectionActions from "src/views/Collection/components/CollectionActions";

const Collection: FC = () => {
  const isCurator = get(FEATURES.CURATOR) === BOOLEAN.TRUE;
  const router = useRouter();
  const { params, tombstoned_dataset_id } = router.query;

  const [hasShownWithdrawToast, setHasShownWithdrawToast] = useState(false);
  const [isUploadingLink, setIsUploadingLink] = useState(false);

  let id = "";

  if (Array.isArray(params)) {
    id = params[0];
  } else if (params) {
    id = params;
  }

  const { data: collection, isError, isFetching } = useCollection({ id });

  useEffect(() => {
    if (
      hasShownWithdrawToast ||
      !tombstoned_dataset_id ||
      !collection ||
      isTombstonedCollection(collection)
    )
      return;

    Toast.show({
      icon: IconNames.ISSUE,
      intent: Intent.PRIMARY,
      message:
        "A dataset was withdrawn. You've been redirected to the parent collection.",
    });
    removeParams("tombstoned_dataset_id");
    setHasShownWithdrawToast(true);
  }, [tombstoned_dataset_id, collection, hasShownWithdrawToast]);

  useEffect(() => {
    if (isTombstonedCollection(collection)) {
      const redirectUrl = ROUTES.HOMEPAGE;
      router.push(redirectUrl + "?tombstoned_collection_id=" + id);
    }
  }, [collection, id, router]);

  /* Pop toast if user has come from Explorer with work in progress */
  useExplainNewTab(
    "To maintain your in-progress work on the previous dataset, we opened this collection in a new tab."
  );

  if (!collection || isError || isTombstonedCollection(collection)) {
    return null;
  }

  const hasRevision = isCollectionHasPrivateRevision(collection);
  const isRevision = isCollectionPrivateRevision(collection);

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

  const collectionConsortia = collection.consortia;
  const collectionMetadataLinks = buildCollectionMetadataLinks(
    collection.links,
    collection.summaryCitation,
    collection.contact_name,
    collection.contact_email
  );

  return (
    <>
      <Head>
        <title>{collection.name} - CZ CELLxGENE Discover</title>
      </Head>
      <CollectionView>
        {/* Collection revision status callout */}
        <CollectionRevisionStatusCallout collection={collection} />
        {/* Collection title and actions */}
        <CollectionHero>
          <h3 data-testid="collection-name">{collection.name}</h3>
          {/* Collection actions */}
          <CollectionActions
            collection={collection}
            hasRevision={hasRevision}
            isPublishable={isPublishable}
            isRevision={isRevision}
            setIsUploadingLink={setIsUploadingLink}
          />
        </CollectionHero>
        {/* Collection consortia, description and metadata */}
        <CollectionDetail>
          {collectionConsortia.length > 0 && (
            <CollectionConsortia>
              {collectionConsortia.join(", ")}
            </CollectionConsortia>
          )}
          <CollectionDescription description={collection.description} />
          {collectionMetadataLinks.length > 0 && (
            <CollectionMetadata
              collectionMetadataLinks={collectionMetadataLinks}
            />
          )}
        </CollectionDetail>
        {/* Collection grid */}
        {/* TODO Reusing DatasetTab as-is as functionality is too dense to refactor for this iteration of filter. Complete refactor (including update to React Table) can be done when filter is productionalized. */}
        <DatasetTab
          collectionId={id}
          datasets={datasets}
          isRevision={isRevision}
          visibility={collection.visibility}
        />
      </CollectionView>
      {/* May be added in the future after sign off */}
      {/* <BottomBanner /> */}
    </>
  );
};

export default Collection;
