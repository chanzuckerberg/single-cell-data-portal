import { Position } from "@blueprintjs/core";
import { useMemo } from "react";
import { Collection } from "src/common/entities";
import RawMoreDropdown from "src/components/common/MoreDropdown";
import Menu from "./components/Menu";
import { DeleteCollectionFn } from "src/views/Collection/components/CollectionActions";

interface Props {
  collection: Collection;
  handleDeleteCollection: DeleteCollectionFn;
  isDeleting: boolean;
  isReorderUX: boolean;
  isRevision: boolean;
}

const MoreDropdown = ({
  collection,
  handleDeleteCollection,
  isDeleting,
  isReorderUX,
  isRevision,
}: Props) => {
  const popoverProps = useMemo(() => {
    return {
      content: (
        <Menu
          collection={collection}
          handleDeleteCollection={handleDeleteCollection}
          isDeleting={isDeleting}
          isReorderUX={isReorderUX}
          isRevision={isRevision}
        />
      ),
      position: Position.BOTTOM,
    };
  }, [collection, handleDeleteCollection, isDeleting, isReorderUX, isRevision]);

  return <RawMoreDropdown popoverProps={popoverProps} />;
};

export default MoreDropdown;
