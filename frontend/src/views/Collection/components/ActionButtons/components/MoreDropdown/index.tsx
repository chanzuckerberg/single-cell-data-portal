import { Position } from "@blueprintjs/core";
import React, { useMemo } from "react";
import RawMoreDropdown from "src/components/common/MoreDropdown";
import Menu from "./components/Menu";

interface Props {
  id: string;
}

const MoreDropdown = ({ id = "" }: Props) => {
  const popoverProps = useMemo(() => {
    return {
      content: <Menu id={id} />,
      position: Position.BOTTOM,
    };
  }, [id]);

  return <RawMoreDropdown popoverProps={popoverProps} />;
};

export default MoreDropdown;
