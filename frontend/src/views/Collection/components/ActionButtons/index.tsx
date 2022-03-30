import { Collection } from "src/common/entities";
import { FEATURES } from "src/common/featureFlags/features";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
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
  const isFilterEnabled = useFeatureFlag(FEATURES.FILTER);
  const Actions = isFilterEnabled ? CollectionActions : Wrapper;
  return (
    <Actions>
      <MoreDropdown id={id} isRevision={!!revisionOf} visibility={visibility} />
      <AddButton addNewFile={addNewFile} />
      <PublishCollection
        id={id}
        isPublishable={isPublishable}
        revisionOf={revisionOf}
      />
    </Actions>
  );
};

export default ActionButtons;
