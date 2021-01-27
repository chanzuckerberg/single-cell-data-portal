import { Button as RawButton, H6, Intent } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import React, { FC, useState } from "react";
import { useDeleteDataset } from "src/common/queries/datasets";
import Alert from "src/components/Alert";

interface Props {
  id?: string;
  collectionId?: string;
  Button?: React.ElementType;
}

const DeleteDataset: FC<Props> = ({
  id,
  collectionId,
  Button = DefaultButton,
}) => {
  const [isOpen, setIsOpen] = useState(false);

  const [deleteDataset] = useDeleteDataset(collectionId);

  const toggleAlert = () => {
    setIsOpen(!isOpen);
  };

  return (
    <>
      <Button disabled={!id} onClick={toggleAlert} />
      <Alert
        cancelButtonText={"Cancel"}
        confirmButtonText={"Delete Dataset"}
        intent={Intent.DANGER}
        isOpen={isOpen}
        onCancel={toggleAlert}
        onConfirm={() => {
          deleteDataset(id);
          toggleAlert();
        }}
      >
        <H6>Are you sure you want to delete this dataset?</H6>
        <p>You cannot undo this action</p>
      </Alert>
    </>
  );
};

function DefaultButton({ ...props }) {
  return <RawButton icon={IconNames.TRASH} minimal {...props} />;
}

export default DeleteDataset;
