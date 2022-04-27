import { Intent, Popover, Position } from "@blueprintjs/core";
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
}

const AddButton = ({ addNewFile }: Props) => {
  return (
    <Popover
      position={Position.BOTTOM_LEFT}
      content={<Content addNewFile={addNewFile} />}
    >
      <StyledOutlineButton intent={Intent.PRIMARY} outlined>
        Add
      </StyledOutlineButton>
    </Popover>
  );
};

export default AddButton;
