import { Position } from "@blueprintjs/core";
import { useMemo } from "react";
import { Collection } from "src/common/entities";
import RawMoreDropdown from "src/components/common/MoreDropdown";
import Menu from "./components/Menu";

interface Props {
  collection: Collection;
  handleDeleteCollection: () => void;
  isDeleting: boolean;
  isRevision: boolean;
}

const MoreDropdown = ({
  collection,
  handleDeleteCollection,
  isDeleting,
  isRevision,
}: Props) => {
  const popoverProps = useMemo(() => {
    return {
      content: (
        <Menu
          collection={collection}
          handleDeleteCollection={handleDeleteCollection}
          isDeleting={isDeleting}
          isRevision={isRevision}
        />
      ),
      position: Position.BOTTOM,
    };
  }, [collection, handleDeleteCollection, isDeleting, isRevision]);

  return <RawMoreDropdown popoverProps={popoverProps} />;
};

export default MoreDropdown;
