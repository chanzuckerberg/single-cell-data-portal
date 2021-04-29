import {
  IMenuItemProps,
  Intent,
  Menu as RawMenu,
  MenuItem,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import React from "react";
import DropboxChooser, {
  Props as ChooserProps,
} from "src/components/DropboxChooser";
import DeleteDataset from "../../../DeleteDataset";

const DeleteButton = (props: IMenuItemProps) => {
  return (
    <MenuItem
      {...props}
      shouldDismissPopover={false}
      icon={IconNames.TRASH}
      intent={Intent.DANGER}
      text="Delete Dataset"
    />
  );
};

const UpdateButton = (props: Partial<IMenuItemProps>) => {
  return (
    <MenuItem
      {...props}
      shouldDismissPopover={false}
      icon={IconNames.EDIT}
      intent={Intent.NONE}
      text="Update Dataset"
    />
  );
};

interface Props {
  collectionId?: string;
  datasetId?: string;
  isRevision: boolean;
  onUploadFile: ChooserProps["onUploadFile"];
}

const Menu = ({
  collectionId = "",
  datasetId = "",
  isRevision,
  onUploadFile,
}: Props) => {
  return (
    <RawMenu>
      <DeleteDataset
        collectionId={collectionId}
        id={datasetId}
        Button={DeleteButton}
      />
      {isRevision && (
        <DropboxChooser onUploadFile={onUploadFile}>
          <UpdateButton />
        </DropboxChooser>
      )}
    </RawMenu>
  );
};

export default Menu;
