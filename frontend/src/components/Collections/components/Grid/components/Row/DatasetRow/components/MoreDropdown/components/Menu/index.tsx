import {
  IMenuItemProps,
  Intent,
  Menu as RawMenu,
  MenuItem,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import DropboxChooser, {
  Props as ChooserProps,
} from "src/components/DropboxChooser";
import styled from "styled-components";
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
  revisionsEnabled: boolean;
  onUploadFile: ChooserProps["onUploadFile"];
  isLoading: boolean;
}

const StyledMenu = styled(RawMenu)`
  border-radius: 3px;
  padding: 8px;
  & > li:last-child {
    margin-bottom: 0;
  }
`;

const Menu = ({
  collectionId = "",
  datasetId = "",
  revisionsEnabled,
  onUploadFile,
  isLoading,
}: Props): JSX.Element => {
  return (
    <StyledMenu>
      {revisionsEnabled && (
        <DropboxChooser onUploadFile={onUploadFile}>
          <UpdateButton disabled={isLoading} />
        </DropboxChooser>
      )}
      <DeleteDataset
        collectionId={collectionId}
        id={datasetId}
        Button={DeleteButton}
      />
    </StyledMenu>
  );
};

export default Menu;
