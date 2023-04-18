import { ACCESS_TYPE, Collection, VISIBILITY_TYPE } from "src/common/entities";
import PublishCollection from "src/components/Collections/components/PublishCollection";
import { Props as DropboxChooserProps } from "src/components/DropboxChooser";
import AddButton from "./components/AddButton";
import MoreDropdown from "./components/MoreDropdown";
import { CollectionActions as Actions } from "./style";
import DeleteCollectionButton from "src/views/Collection/components/CollectionActions/components/DeleteButton";
import CreateRevisionButton from "src/views/Collection/components/CollectionActions/components/CreateRevisionButton";
import React, { useState } from "react";
import {
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

interface Props {
  addNewFile: DropboxChooserProps["onUploadFile"];
  collection: Collection;
  handleDeleteCollection: () => void;
  hasRevision: boolean;
  isDeleting: boolean;
  isPublishable: boolean;
  isRevision: boolean;
}

const CollectionActions = ({
  addNewFile,
  collection,
  handleDeleteCollection,
  hasRevision,
  isDeleting,
  isPublishable,
  isRevision,
}: Props): JSX.Element | null => {
  const { id, name, revision_of } = collection;
  const [isPublishOpen, setIsPublishOpen] = useState(false);
  const router = useRouter();
  const { mutateAsync: createRevision } = useCreateRevision();
  const { mutateAsync: deleteCollection, isLoading: isDeletingCollection } =
    useDeleteCollection();
  const { mutateAsync: publishCollection, isLoading: isPublishing } =
    usePublishCollection();

  // Creates a revision of the collection and routes to the private revision collection.
  const handleCreateRevision = async (): Promise<void> => {
    await createRevision(id, {
      onSuccess: (revision) => {
        router.push(ROUTES.COLLECTION.replace(":id", revision.id));
      },
    });
  };

  // Deletes a private collection and routes to the collections index.
  const handleDeletePrivateCollection = async (): Promise<void> => {
    await deleteCollection(collection, {
      onSuccess: () => {
        console.log("Successfully deleted private collection!");
        router.push(ROUTES.COLLECTIONS);
      },
    });
  };

  // Deletes a private revision and routes to the published collection.
  const handleDeletePrivateRevisionCollection = async (): Promise<void> => {
    await deleteCollection(collection, {
      onSuccess: () => {
        console.log("Successfully deleted private revision collection!");
        if (revision_of) {
          router.push(ROUTES.COLLECTION.replace(":id", revision_of));
        }
      },
    });
  };

  // Publishes a private collection or private revision and routes the published private revision to the published collection.
  const handlePublishCollection = async () => {
    const payload = JSON.stringify({
      data_submission_policy_version: POLICY_BULLETS.version,
    });
    await publishCollection(
      { collection, payload },
      {
        onSuccess: () => {
          if (revision_of) {
            Toast.show({
              icon: IconNames.TICK,
              intent: Intent.SUCCESS,
              message: "New version published",
            });
            router.push(ROUTES.COLLECTION.replace(":id", revision_of));
          }
        },
      }
    );
    setIsPublishOpen(false);
  };

  return collection.access_type === ACCESS_TYPE.WRITE ? (
    <Actions data-testid="collection-actions">
      {/* Collection is either private, or a private revision */}
      {collection.visibility === VISIBILITY_TYPE.PRIVATE && (
        <>
          <MoreDropdown
            collection={collection}
            handleDeleteCollection={
              isRevision
                ? handleDeletePrivateRevisionCollection
                : handleDeletePrivateCollection
            }
            isDeleting={isDeletingCollection}
            isRevision={isRevision}
          />
          <AddButton addNewFile={addNewFile} />
          <PublishCollection
            handlePublishCollection={handlePublishCollection}
            isPublishable={isPublishable}
            isPublishing={isPublishing}
            isPublishOpen={isPublishOpen}
            revision_of={revision_of}
            setIsPublishOpen={setIsPublishOpen}
          />
        </>
      )}
      {/* Collection is public */}
      {collection.visibility === VISIBILITY_TYPE.PUBLIC && (
        <>
          <DeleteCollectionButton
            disabled={hasRevision}
            collectionName={name}
            handleConfirm={handleDeleteCollection}
            loading={isDeleting}
          />
          {!hasRevision && (
            <CreateRevisionButton handleCreateRevision={handleCreateRevision} />
          )}
        </>
      )}
    </Actions>
  ) : null;
};

export default CollectionActions;
