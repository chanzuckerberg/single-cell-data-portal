import { Position } from "@blueprintjs/core";
import React, { useMemo } from "react";
import RawMoreDropdown from "src/components/common/MoreDropdown";
import Menu from "./components/Menu";

interface Props {
  collectionId?: string;
  datasetId?: string;
  isRevision: boolean;
}

const MoreDropdown = ({
  collectionId = "",
  datasetId = "",
  isRevision,
}: Props) => {
  const popoverProps = useMemo(() => {
    return {
      content: (
        <Menu
          collectionId={collectionId}
          datasetId={datasetId}
          isRevision={isRevision}
        />
      ),
      position: Position.BOTTOM,
    };
  }, [collectionId, datasetId, isRevision]);

  return <RawMoreDropdown popoverProps={popoverProps} />;
};

export default MoreDropdown;
