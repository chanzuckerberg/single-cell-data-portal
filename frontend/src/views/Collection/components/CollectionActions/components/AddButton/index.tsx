import { Popover, Position } from "@blueprintjs/core";
import {
  Props as DropboxChooserProps,
  UploadingFile,
} from "src/components/DropboxChooser";
import Content from "./components/Content";
import { ActionButton as Button } from "src/views/Collection/components/CollectionActions/style";

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
      <Button sdsStyle="square" sdsType="secondary">
        Add
      </Button>
    </Popover>
  );
};

export default AddButton;
