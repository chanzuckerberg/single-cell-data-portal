import { Button as RawButton, H6, Intent } from "@blueprintjs/core";
import loadable from "@loadable/component";
import { useRouter } from "next/router";
import * as React from "react";
import { FC, useState } from "react";
import { ROUTES } from "src/common/constants/routes";
import { Collection } from "src/common/entities";
import {
  useDeleteCollection,
  useDeletePrivateRevisionCollection,
} from "src/common/queries/collections";

const AsyncAlert = loadable(
  () =>
    /*webpackChunkName: 'src/components/Alert' */ import("src/components/Alert")
);

interface Props {
  Button?: React.ElementType;
  collection: Collection;
  isRevision: boolean;
}

const DeleteCollection: FC<Props> = ({
  Button = RawButton,
  collection,
  isRevision,
}) => {
  const { id, revision_of } = collection;
  const router = useRouter();
  const {
    mutateAsync: deletePrivateCollection,
    isLoading: isDeletePrivateCollectionLoading,
  } = useDeleteCollection(id);
  const {
    mutateAsync: deletePrivateRevision,
    isLoading: isDeletePrivateRevisionLoading,
  } = useDeletePrivateRevisionCollection(collection);

  const handleDeletePrivateCollection = async (): Promise<void> => {
    // (thuang): `deleteMutation` should have supported `onSuccess` callback,
    // but it doesn't seem to work. Thus we await the promise then redirect!
    await deletePrivateCollection(id);

    router.push(ROUTES.COLLECTIONS);
  };

  // Deletes private revision collection (cancel the revision) and routes to the published collection.
  const handleDeletePrivateRevision = async (): Promise<void> => {
    // (thuang): `deleteMutation` should have supported `onSuccess` callback,
    // but it doesn't seem to work. Thus we await the promise then redirect!
    await deletePrivateRevision(id);

    if (revision_of) {
      router.push(ROUTES.COLLECTION.replace(":id", revision_of));
    }
  };

  const [isOpen, setIsOpen] = useState(false);

  const toggleAlert = () => {
    setIsOpen(!isOpen);
  };

  const handleHover = () => {
    AsyncAlert.preload();
  };

  return (
    <>
      <Button
        onClick={toggleAlert}
        onMouseEnter={handleHover}
        isRevision={isRevision}
      />

      {isOpen && (
        <AsyncAlert
          cancelButtonText={"Cancel"}
          confirmButtonText={
            isRevision ? "Cancel Revision" : "Delete Collection"
          }
          intent={Intent.DANGER}
          isOpen={isOpen}
          onCancel={toggleAlert}
          onConfirm={
            isRevision
              ? handleDeletePrivateRevision
              : handleDeletePrivateCollection
          }
          loading={
            isDeletePrivateCollectionLoading || isDeletePrivateRevisionLoading
          }
        >
          <H6>
            Are you sure you want to {isRevision ? "cancel" : "delete"} this{" "}
            {isRevision ? "revision" : "collection"}?
          </H6>
          <p>
            {isRevision
              ? "Cancelling this revision"
              : "Deleting this collection"}{" "}
            will also delete any uploaded datasets.{" "}
            {!isRevision &&
              "If you’ve shared this collection or its datasets with anyone, they will also lose access. "}
            You cannot undo this action.
          </p>
        </AsyncAlert>
      )}
    </>
  );
};

export default DeleteCollection;
