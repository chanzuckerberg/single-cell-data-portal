import { Position } from "@blueprintjs/core";
import { useMemo } from "react";
import { Collection } from "src/common/entities";
import RawMoreDropdown from "src/components/common/MoreDropdown";
import Menu from "./components/Menu";

interface Props {
  id: string;
  isRevision: boolean;
  visibility: Collection["visibility"];
}

const MoreDropdown = ({ id = "", isRevision, visibility }: Props) => {
  const popoverProps = useMemo(() => {
    return {
      content: <Menu id={id} isRevision={isRevision} visibility={visibility} />,
      position: Position.BOTTOM,
    };
  }, [id, isRevision]);

  return <RawMoreDropdown popoverProps={popoverProps} />;
};

export default MoreDropdown;
