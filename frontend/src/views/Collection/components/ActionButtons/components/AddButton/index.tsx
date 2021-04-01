import { Button, Intent, Popover, Position } from "@blueprintjs/core";
import React from "react";
import {
  Props as DropboxChooserProps,
  UploadingFile,
} from "src/components/DropboxChooser";
import Content from "./components/Content";

export interface UploadedFiles {
  [datasetID: string]: UploadingFile;
}

interface Props {
  addNewFile: DropboxChooserProps["onUploadFile"];
}

const AddButton = ({ addNewFile }: Props) => {
  return (
    <Popover
      position={Position.BOTTOM_LEFT}
      content={<Content addNewFile={addNewFile} />}
    >
      <Button intent={Intent.PRIMARY} outlined>
        Add
      </Button>
    </Popover>
  );
};

export default AddButton;
