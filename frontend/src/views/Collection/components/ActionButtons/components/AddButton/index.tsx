import { Button, Intent, Popover, Position } from "@blueprintjs/core";
import { StyledOutlineButton } from "src/components/common/Button/common/style";
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
  isFilterEnabled?: boolean;
}

const AddButton = ({ addNewFile, isFilterEnabled }: Props) => {
  const AddDatasetButton = isFilterEnabled ? StyledOutlineButton : Button;
  return (
    <Popover
      position={Position.BOTTOM_LEFT}
      content={<Content addNewFile={addNewFile} />}
    >
      <AddDatasetButton intent={Intent.PRIMARY} outlined>
        Add
      </AddDatasetButton>
    </Popover>
  );
};

export default AddButton;
