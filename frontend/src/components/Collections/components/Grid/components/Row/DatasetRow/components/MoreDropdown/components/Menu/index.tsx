import {
  IMenuItemProps,
  Intent,
  Menu as RawMenu,
  MenuItem,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import React from "react";
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

interface Props {
  collectionId?: string;
  datasetId?: string;
}

const Menu = ({ collectionId = "", datasetId = "" }: Props) => {
  return (
    <RawMenu>
      <DeleteDataset
        collectionId={collectionId}
        id={datasetId}
        Button={DeleteButton}
      />
    </RawMenu>
  );
};

export default Menu;
