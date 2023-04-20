import { Collection } from "src/common/entities";
import PublishCollection from "src/components/Collections/components/PublishCollection";
import {
  Props as DropboxChooserProps,
  UploadingFile,
} from "src/components/DropboxChooser";
import AddButton from "./components/AddButton";
import MoreDropdown from "./components/MoreDropdown";
import { CollectionActions } from "./style";
export interface UploadedFiles {
  [datasetID: string]: UploadingFile;
}

interface Props {
  addNewFile: DropboxChooserProps["onUploadFile"];
  id: Collection["id"];
  isPublishable: boolean;
  revisionOf: Collection["revision_of"];
  visibility: Collection["visibility"];
}

const ActionButtons = ({
  addNewFile,
  id,
  isPublishable,
  revisionOf,
  visibility,
}: Props): JSX.Element => {
  return (
    <CollectionActions data-testid="collection-actions">
      <MoreDropdown id={id} isRevision={!!revisionOf} visibility={visibility} />
      <AddButton addNewFile={addNewFile} />
      <PublishCollection
        id={id}
        isPublishable={isPublishable}
        revisionOf={revisionOf}
      />
    </CollectionActions>
  );
};

export default ActionButtons;
