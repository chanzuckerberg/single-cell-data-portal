import React from "react";
import { Collection } from "src/common/entities";
import PublishCollection from "src/components/Collections/components/PublishCollection";
import {
  Props as DropboxChooserProps,
  UploadingFile,
} from "src/components/DropboxChooser";
import AddButton from "./components/AddButton";
import MoreDropdown from "./components/MoreDropdown";
import { Wrapper } from "./style";
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
  return (
    <Wrapper>
      <MoreDropdown id={id} isRevision={isRevision} />

      <AddButton addNewFile={addNewFile} />
      <PublishCollection
        id={id}
        isPublishable={isPublishable}
        isRevision={isRevision}
      />
    </Wrapper>
  );
};

export default ActionButtons;
