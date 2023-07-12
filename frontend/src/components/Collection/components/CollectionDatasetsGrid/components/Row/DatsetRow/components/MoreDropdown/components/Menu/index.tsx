import {
  IMenuItemProps,
  Intent,
  Menu as RawMenu,
  MenuItem,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import styled from "@emotion/styled";
import DropboxChooser, {
  Props as ChooserProps,
} from "src/components/DropboxChooser";
import DeleteDataset from "../../../DeleteDataset";
import { Collection } from "src/common/entities";

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
  collectionId: Collection["id"];
  datasetId?: string;
  isPublished: boolean;
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
  collectionId,
  datasetId = "",
  isPublished,
  revisionsEnabled,
  onUploadFile,
  isLoading,
}: Props): JSX.Element => {
  // A dataset may be deleted if the collection is private, or the dataset has not been previously
  // published; where the published_at property is used to determine whether the dataset has been previously published.
  const shouldShowDelete = !revisionsEnabled || !isPublished;
  return (
    <StyledMenu>
      {revisionsEnabled && (
        <DropboxChooser onUploadFile={onUploadFile}>
          <UpdateButton disabled={isLoading} />
        </DropboxChooser>
      )}
      {shouldShowDelete && (
        <DeleteDataset
          Button={DeleteButton}
          collectionId={collectionId}
          datasetId={datasetId}
        />
      )}
    </StyledMenu>
  );
};

export default Menu;
