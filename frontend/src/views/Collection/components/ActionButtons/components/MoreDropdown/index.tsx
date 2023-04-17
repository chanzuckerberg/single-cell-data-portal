import { Position } from "@blueprintjs/core";
import { useMemo } from "react";
import { Collection } from "src/common/entities";
import RawMoreDropdown from "src/components/common/MoreDropdown";
import Menu from "./components/Menu";

interface Props {
  collection: Collection;
  isRevision: boolean;
}

const MoreDropdown = ({ collection, isRevision }: Props) => {
  const popoverProps = useMemo(() => {
    return {
      content: <Menu collection={collection} isRevision={isRevision} />,
      position: Position.BOTTOM,
    };
  }, [collection, isRevision]);

  return <RawMoreDropdown popoverProps={popoverProps} />;
};

export default MoreDropdown;
