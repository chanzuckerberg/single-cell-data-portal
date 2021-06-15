import { Intent } from "@blueprintjs/core";
import React, { ReactElement } from "react";
import { Dataset } from "src/common/entities";
import { StyledTag } from "./style";

interface Props {
  dataset: Dataset;
}

enum REVISION_STATUS {
  DELETED = "Deleted",
  UPDATED = "Updated",
  NEW = "New",
}

const statusToIntent = {
  [REVISION_STATUS.DELETED]: Intent.DANGER,
  [REVISION_STATUS.UPDATED]: Intent.PRIMARY,
  [REVISION_STATUS.NEW]: Intent.SUCCESS,
};

const RevisionStatusTag = ({ dataset }: Props): ReactElement | null => {
  let revisionStatus = REVISION_STATUS.NEW;
  const isPublished = dataset.published;
  if (!isPublished) {
    if (dataset.original_id) {
      revisionStatus = REVISION_STATUS.UPDATED;
    } else if (!dataset.original_id) {
      revisionStatus = REVISION_STATUS.NEW;
    }
  } else if (dataset.tombstone) {
    revisionStatus = REVISION_STATUS.DELETED;
  } else {
    return null;
  }
  return (
    <StyledTag minimal intent={statusToIntent[revisionStatus]}>
      {revisionStatus}
    </StyledTag>
  );
};

export default RevisionStatusTag;
