import { Collection } from "src/common/entities";
import PublishCollection from "src/components/Collections/components/PublishCollection";
import {
  Props as DropboxChooserProps,
  UploadingFile,
} from "src/components/DropboxChooser";
import AddButton from "./components/AddButton";
import MoreDropdown from "./components/MoreDropdown";
import { CollectionActions, Wrapper } from "./style";
export interface UploadedFiles {
  [datasetID: string]: UploadingFile;
}

interface Props {
  addNewFile: DropboxChooserProps["onUploadFile"];
  id: Collection["id"];
  isFilterEnabled?: boolean;
  isPublishable: boolean;
  isRevision: boolean;
}

const ActionButtons = ({
  addNewFile,
  id,
  isFilterEnabled = false,
  isPublishable,
  isRevision,
}: Props): JSX.Element => {
  const Actions = isFilterEnabled ? CollectionActions : Wrapper;
  return (
    <Actions>
      <MoreDropdown id={id} isRevision={isRevision} />
      <AddButton addNewFile={addNewFile} isFilterEnabled={isFilterEnabled} />
      <PublishCollection
        id={id}
        isFilterEnabled={isFilterEnabled}
        isPublishable={isPublishable}
        isRevision={isRevision}
      />
    </Actions>
  );
};

export default ActionButtons;
