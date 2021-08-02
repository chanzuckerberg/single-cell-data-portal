import { Button as RawButton, H6, Intent } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import loadable from "@loadable/component";
import React, { FC, useState } from "react";
import { useDeleteDataset } from "src/common/queries/datasets";

const AsyncAlert = loadable(
  () =>
    /*webpackChunkName: 'src/components/Alert' */ import("src/components/Alert")
);

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
  const [isLoading, setIsLoading] = useState(false);

  const [deleteDataset] = useDeleteDataset(collectionId);

  const toggleAlert = () => {
    setIsOpen(!isOpen);
  };

  const handleHover = () => {
    AsyncAlert.preload();
  };

  return (
    <>
      <Button onMouseEnter={handleHover} disabled={!id} onClick={toggleAlert} />
      {isOpen && (
        <AsyncAlert
          loading={isLoading}
          cancelButtonText="Cancel"
          confirmButtonText="Delete Dataset"
          intent={Intent.DANGER}
          isOpen={isOpen}
          onCancel={toggleAlert}
          onConfirm={async () => {
            setIsLoading(true);
            await deleteDataset(id);
            toggleAlert();
          }}
        >
          <H6>Are you sure you want to delete this dataset?</H6>
          <p>You cannot undo this action</p>
        </AsyncAlert>
      )}
    </>
  );
};

function DefaultButton({ ...props }) {
  return <RawButton icon={IconNames.TRASH} minimal {...props} />;
}

export default DeleteDataset;
