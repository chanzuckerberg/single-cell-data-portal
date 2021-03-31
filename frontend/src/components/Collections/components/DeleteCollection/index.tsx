import { Button, H6, Intent } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import loadable from "@loadable/component";
import { useRouter } from "next/router";
import React, { FC, useState } from "react";
import { ROUTES } from "src/common/constants/routes";
import { Collection } from "src/common/entities";
import { useDeleteCollection } from "src/common/queries/collections";

const AsyncAlert = loadable(
  () =>
    /*webpackChunkName: 'src/components/Alert' */ import("src/components/Alert")
);
interface Props {
  id: Collection["id"];
}

const DeleteCollection: FC<Props> = ({ id }) => {
  const [deleteMutation, { isSuccess }] = useDeleteCollection(id);
  const router = useRouter();

  if (isSuccess) router.push(ROUTES.MY_COLLECTIONS);

  const handleDelete = () => {
    deleteMutation(id);
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
        intent={Intent.DANGER}
        minimal
        text="Delete"
        icon={IconNames.TRASH}
        onClick={toggleAlert}
        onMouseEnter={handleHover}
      />

      {isOpen && (
        <AsyncAlert
          cancelButtonText={"Cancel"}
          confirmButtonText={"Delete Collection"}
          intent={Intent.DANGER}
          isOpen={isOpen}
          onCancel={toggleAlert}
          onConfirm={() => {
            handleDelete();
            toggleAlert();
          }}
        >
          <H6>Are you sure you want to delete this collection?</H6>
          <p>
            Deleting this collection will also delete any uploaded datasets. If
            you’ve shared this collection or its datasets with anyone, they will
            also lose access. You cannot undo this action.
          </p>
        </AsyncAlert>
      )}
    </>
  );
};

export default DeleteCollection;
