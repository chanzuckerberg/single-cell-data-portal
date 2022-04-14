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
  isRevision: boolean;
}

const ActionButtons = ({
  addNewFile,
  id,
  isPublishable,
  isRevision,
}: Props): JSX.Element => {
  const isFilterEnabled = useFeatureFlag(FEATURES.FILTER);
  const Actions = isFilterEnabled ? CollectionActions : Wrapper;
  return (
    <Actions>
      <MoreDropdown id={id} isRevision={isRevision} />
      <AddButton addNewFile={addNewFile} />
      <PublishCollection
        id={id}
        isPublishable={isPublishable}
        isRevision={isRevision}
      />
    </Actions>
  );
};

export default ActionButtons;
