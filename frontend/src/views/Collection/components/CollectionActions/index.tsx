import { Collection, VISIBILITY_TYPE } from "src/common/entities";
import PublishCollection from "src/components/Collections/components/PublishCollection";
import {
  Props as DropboxChooserProps,
  UploadingFile,
} from "src/components/DropboxChooser";
import AddButton from "./components/AddButton";
import MoreDropdown from "./components/MoreDropdown";
import { CollectionActions as Actions } from "./style";
import DeleteCollectionButton from "src/views/Collection/components/CollectionActions/components/DeleteButton";
import CreateRevisionButton from "src/views/Collection/components/CollectionActions/components/CreateRevisionButton";
import React from "react";
import { useCreateRevision } from "src/common/queries/collections";
import { ROUTES } from "src/common/constants/routes";
import { useRouter } from "next/router";

export interface UploadedFiles {
  [datasetID: string]: UploadingFile;
}

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
}: Props): JSX.Element => {
  const { id, name, revision_of, visibility } = collection;
  const router = useRouter();
  const { mutateAsync: createRevision } = useCreateRevision();

  // Creates a revision of the collection and routes to the private revision collection.
  const handleCreateRevision = async (): Promise<void> => {
    await createRevision(id, {
      onSuccess: (revision) => {
        router.push(ROUTES.COLLECTION.replace(":id", revision.id));
      },
    });
  };

  return (
    <Actions data-testid="collection-actions">
      {/* Collection is either private, or a private revision */}
      {collection.visibility === VISIBILITY_TYPE.PRIVATE && (
        <>
          <MoreDropdown
            id={id}
            isRevision={isRevision}
            visibility={visibility}
          />
          <AddButton addNewFile={addNewFile} />
          <PublishCollection
            id={id}
            isPublishable={isPublishable}
            revisionOf={revision_of}
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
  );
};

export default CollectionActions;
