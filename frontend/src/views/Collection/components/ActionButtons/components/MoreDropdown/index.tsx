import { Position } from "@blueprintjs/core";
import React, { useMemo } from "react";
import RawMoreDropdown from "src/components/common/MoreDropdown";
import Menu from "./components/Menu";

interface Props {
  id: string;
  isRevision: boolean;
}

const MoreDropdown = ({ id = "", isRevision }: Props) => {
  const popoverProps = useMemo(() => {
    return {
      content: <Menu id={id} isRevision={isRevision} />,
      position: Position.BOTTOM,
    };
  }, [id, isRevision]);

  return <RawMoreDropdown popoverProps={popoverProps} />;
};

export default MoreDropdown;
