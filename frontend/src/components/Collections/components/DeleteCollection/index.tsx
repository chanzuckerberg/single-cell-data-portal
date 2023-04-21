import { Button as RawButton, H6, Intent } from "@blueprintjs/core";
import loadable from "@loadable/component";
import * as React from "react";
import { FC, useEffect, useState } from "react";
import { DeleteCollectionFn } from "src/views/Collection/components/CollectionActions";

const AsyncAlert = loadable(
  () =>
    /*webpackChunkName: 'src/components/Alert' */ import("src/components/Alert")
);

interface Props {
  Button?: React.ElementType;
  handleDeleteCollection: DeleteCollectionFn;
  isDeleting: boolean;
  isRevision: boolean;
}

const DeleteCollection: FC<Props> = ({
  Button = RawButton,
  handleDeleteCollection,
  isDeleting,
  isRevision,
}) => {
  const [isOpen, setIsOpen] = useState(false);

  const handleHover = () => {
    AsyncAlert.preload();
  };

  // Closes delete collection dialog when component unmounts.
  // If a private revision collection is successfully deleted or a private collection is deleted, the user will be
  // directed to the corresponding published collection or the collections index, respectively.
  // In either case, the component will unmount and the delete collection dialog will be closed.
  useEffect(() => {
    return () => {
      setIsOpen(false);
    };
  }, []);

  return (
    <>
      <Button
        onClick={() => setIsOpen(true)}
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
          onCancel={() => setIsOpen(false)}
          onConfirm={handleDeleteCollection}
          loading={isDeleting}
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
