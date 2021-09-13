import { Button as RawButton, H6, Intent } from "@blueprintjs/core";
import loadable from "@loadable/component";
import { useRouter } from "next/router";
import * as React from "react";
import { FC, useState } from "react";
import { ROUTES } from "src/common/constants/routes";
import { Collection } from "src/common/entities";
import { useDeleteCollection } from "src/common/queries/collections";

const AsyncAlert = loadable(
  () =>
    /*webpackChunkName: 'src/components/Alert' */ import("src/components/Alert")
);
interface Props {
  id: Collection["id"];
  Button?: React.ElementType;
  isRevision: boolean;
}

const DeleteCollection: FC<Props> = ({
  id,
  Button = RawButton,
  isRevision,
}) => {
  const [deleteMutation] = useDeleteCollection(id);
  const [isLoading, setIsLoading] = useState(false);
  const router = useRouter();

  const handleDelete = async () => {
    setIsLoading(true);
    await deleteMutation(id);
    router.push(ROUTES.MY_COLLECTIONS);
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
          onConfirm={handleDelete}
          loading={isLoading}
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
