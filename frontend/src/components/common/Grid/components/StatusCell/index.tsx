import React, { ElementType, ReactNode } from "react";
import { COLLECTION_STATUS } from "src/common/entities";
import {
  PrivateTag,
  PublishedTag,
  RevisionTag,
  StatusTag,
  StatusTags,
} from "src/components/common/Grid/components/StatusCell/style";

interface Props {
  revisionButton?: ReactNode;
  status?: COLLECTION_STATUS[];
}

/**
 * Returns the status tag element type based on the status.
 * @param status - Collection status.
 * @returns status tag element type.
 */
function getStatusTagElType(status: COLLECTION_STATUS): ElementType {
  switch (status) {
    case COLLECTION_STATUS.REVISION:
      return RevisionTag;
    case COLLECTION_STATUS.PUBLISHED:
      return PublishedTag;
    case COLLECTION_STATUS.PRIVATE:
      return PrivateTag;
  }
}

export default function StatusCell({
  revisionButton,
  status,
}: Props): JSX.Element | null {
  return status ? (
    <StatusTags>
      {status.map((status) => {
        const Tag = getStatusTagElType(status);
        return (
          <StatusTag key={status}>
            <Tag
              color={
                status === COLLECTION_STATUS.PUBLISHED ? "success" : "primary"
              }
              data-testid="status-tag"
              label={status}
              sdsStyle="square"
              sdsType="primary"
            />
            {status === COLLECTION_STATUS.REVISION ? revisionButton : null}
          </StatusTag>
        );
      })}
    </StatusTags>
  ) : null;
}
