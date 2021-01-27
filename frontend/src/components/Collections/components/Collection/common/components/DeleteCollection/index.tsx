import { Button, H6, Intent } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { navigate } from "@reach/router";
import React, { FC, useState } from "react";
import { ROUTES } from "src/common/constants/routes";
import { Collection } from "src/common/entities";
import { useDeleteCollection } from "src/common/queries/collections";
import StyledAlert from "src/components/Collections/common/StyledAlert";

interface Props {
  id: Collection["id"];
}

const DeleteCollection: FC<Props> = ({ id }) => {
  const [deleteMutation] = useDeleteCollection(id);

  const handleDelete = () => {
    deleteMutation(id);
    navigate(ROUTES.MY_COLLECTIONS);
  };

  const [isOpen, setIsOpen] = useState(false);

  const toggleAlert = () => {
    setIsOpen(!isOpen);
  };

  return (
    <>
      <Button
        intent={Intent.DANGER}
        minimal
        text="Delete"
        icon={IconNames.TRASH}
        onClick={toggleAlert}
      />

      <StyledAlert
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
      </StyledAlert>
    </>
  );
};

export default DeleteCollection;
