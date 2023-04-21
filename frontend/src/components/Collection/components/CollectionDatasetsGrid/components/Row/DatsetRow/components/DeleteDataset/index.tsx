import { Button as RawButton, H6, Intent } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import loadable from "@loadable/component";
import * as React from "react";
import { FC, useState } from "react";
import { useDeleteDataset } from "src/common/queries/datasets";
import { Collection } from "src/common/entities";

const AsyncAlert = loadable(
  () =>
    /*webpackChunkName: 'src/components/Alert' */ import("src/components/Alert")
);

interface Props {
  Button?: React.ElementType;
  collectionId: Collection["id"];
  datasetId?: string;
}

const DeleteDataset: FC<Props> = ({
  Button = DefaultButton,
  collectionId,
  datasetId,
}) => {
  const [isOpen, setIsOpen] = useState(false);
  const { mutateAsync: deleteDataset, isLoading } = useDeleteDataset();

  const handleHover = () => {
    AsyncAlert.preload();
  };

  // Deletes dataset.
  const handleDeleteDataset = async (): Promise<void> => {
    if (!datasetId) return;
    await deleteDataset(
      { collectionId, datasetId },
      {
        onSuccess: () => {
          setIsOpen(false);
        },
      }
    );
  };

  return (
    <>
      <Button
        onMouseEnter={handleHover}
        disabled={!datasetId}
        onClick={() => setIsOpen(true)}
      />
      {isOpen && (
        <AsyncAlert
          loading={isLoading}
          cancelButtonText="Cancel"
          confirmButtonText="Delete Dataset"
          intent={Intent.DANGER}
          isOpen={isOpen}
          onCancel={() => setIsOpen(false)}
          onConfirm={handleDeleteDataset}
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
