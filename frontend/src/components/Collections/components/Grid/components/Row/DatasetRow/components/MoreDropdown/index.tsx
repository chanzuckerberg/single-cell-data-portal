import { Position } from "@blueprintjs/core";
import React, { useMemo } from "react";
import RawMoreDropdown from "src/components/common/MoreDropdown";
import Menu from "./components/Menu";

interface Props {
  collectionId?: string;
  datasetId?: string;
}

const MoreDropdown = ({ collectionId = "", datasetId = "" }: Props) => {
  const popoverProps = useMemo(() => {
    return {
      content: <Menu collectionId={collectionId} datasetId={datasetId} />,
      position: Position.BOTTOM,
    };
  }, [collectionId, datasetId]);

  return <RawMoreDropdown popoverProps={popoverProps} />;
};

export default MoreDropdown;
