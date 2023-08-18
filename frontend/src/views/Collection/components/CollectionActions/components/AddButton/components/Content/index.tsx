import { Menu, MenuItem } from "@blueprintjs/core";
import DropboxChooser, {
  Props as DropboxChooserProps,
} from "src/components/DropboxChooser";

interface Props {
  addNewFile: DropboxChooserProps["onUploadFile"];
}

const Content = ({ addNewFile }: Props) => {
  return (
    <Menu>
      <DropboxChooser onUploadFile={addNewFile}>
        <MenuItem text="Add Dataset" />
      </DropboxChooser>
    </Menu>
  );
};

export default Content;
