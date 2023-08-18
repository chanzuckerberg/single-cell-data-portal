import { ACCESS_TYPE, Collection, VISIBILITY_TYPE } from "src/common/entities";
import PublishCollection from "src/components/Collections/components/PublishCollection";
import { UploadingFile } from "src/components/DropboxChooser";
import AddButton from "./components/AddButton";
import MoreDropdown from "./components/MoreDropdown";
import { CollectionActions as Actions } from "./style";
import CreateRevisionButton from "src/views/Collection/components/CollectionActions/components/CreateRevisionButton";
import React, { Dispatch, SetStateAction } from "react";
import {
  useCollectionUploadLinks,
  useCreateRevision,
  useDeleteCollection,
  usePublishCollection,
} from "src/common/queries/collections";
import { ROUTES } from "src/common/constants/routes";
import { useRouter } from "next/router";
import { POLICY_BULLETS } from "src/components/Collections/components/PublishCollection/components/Policy";
import Toast from "src/views/Collection/components/Toast";
import { IconNames } from "@blueprintjs/icons";
import { Intent } from "@blueprintjs/core";

export type CreateRevisionFn = () => void;
export type DeleteCollectionFn = () => void;
export type PublishCollectionFn = () => void;

interface Props {
  collection: Collection;
  hasRevision: boolean;
  isPublishable: boolean;
  isRevision: boolean;
  setIsUploadingLink: Dispatch<SetStateAction<boolean>>;
}

const CollectionActions = ({
  collection,
  hasRevision,
  isPublishable,
  isRevision,
  setIsUploadingLink,
}: Props): JSX.Element | null => {
  const { id, revision_of } = collection;
  const createRevisionMutation = useCreateRevision();
  const uploadLinksMutation = useCollectionUploadLinks(id);
  const deleteCollectionMutation = useDeleteCollection();
  const publishCollectionMutation = usePublishCollection();
  const router = useRouter();

  // Adds a new file to the collection.
  const handleAddNewFile = async (newFile: UploadingFile): Promise<void> => {
    if (!newFile.link) return;

    const payload = JSON.stringify({ url: newFile.link });
    setIsUploadingLink(true);

    await uploadLinksMutation.mutateAsync(
      { collectionId: id, payload },
      {
        onSettled: () => {
          setIsUploadingLink(false);
        },
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

  // Creates a revision of the collection.
  const handleCreateRevision = async (): Promise<void> => {
    await createRevisionMutation.mutateAsync(id, {
      onSuccess: (revision) => {
        // Route to the private revision.
        router.push(ROUTES.COLLECTION.replace(":id", revision.id));
      },
    });
  };

  // Deletes a private, or private revision collection.
  const handleDeleteCollection = async (): Promise<void> => {
    await deleteCollectionMutation.mutateAsync(collection, {
      onSuccess: () => {
        if (revision_of) {
          // If the collection is a private revision, route to the published collection.
          router.push(ROUTES.COLLECTION.replace(":id", revision_of));
        } else {
          // Otherwise, route to the collections page.
          router.push(ROUTES.COLLECTIONS);
        }
      },
    });
  };

  // Publishes a private collection or private revision.
  const handlePublishCollection = async (): Promise<void> => {
    const payload = JSON.stringify({
      data_submission_policy_version: POLICY_BULLETS.version,
    });
    await publishCollectionMutation.mutateAsync(
      { collection, payload },
      {
        onSuccess: () => {
          if (revision_of) {
            // If the collection is a private revision show toast.
            Toast.show({
              icon: IconNames.TICK,
              intent: Intent.SUCCESS,
              message: "New version published",
            });
            // Route to the published collection.
            router.push(ROUTES.COLLECTION.replace(":id", revision_of));
          }
          // Note, there is no need to route if the collection is private.
          // The invalidation of the cache will cause the collection to be re-fetched and
          // collection visibility etc. will be updated.
        },
      }
    );
  };

  return collection.access_type === ACCESS_TYPE.WRITE ? (
    <Actions data-testid="collection-actions">
      {/* Collection is either private, or a private revision */}
      {collection.visibility === VISIBILITY_TYPE.PRIVATE && (
        <>
          <MoreDropdown
            collection={collection}
            handleDeleteCollection={handleDeleteCollection}
            isDeleting={deleteCollectionMutation.isLoading}
            isRevision={isRevision}
          />
          <AddButton addNewFile={handleAddNewFile} />
          <PublishCollection
            handlePublishCollection={handlePublishCollection}
            isPublishable={isPublishable}
            isPublishing={publishCollectionMutation.isLoading}
            revision_of={revision_of}
          />
        </>
      )}
      {/* Collection is public */}
      {collection.visibility === VISIBILITY_TYPE.PUBLIC && (
        <>
          {!hasRevision && (
            <CreateRevisionButton handleCreateRevision={handleCreateRevision} />
          )}
        </>
      )}
    </Actions>
  ) : null;
};

export default CollectionActions;
