import { Fragment, MouseEvent, useState } from "react";
import { Collection } from "src/common/entities";
import Menu from "./components/Menu";
import { DeleteCollectionFn } from "src/views/Collection/components/CollectionActions";
import { Reorder } from "src/views/Collection/hooks/useReorder/common/entities";
import { Button } from "src/views/Collection/components/CollectionActions/components/MoreDropdown/style";

interface Props {
  collection: Collection;
  handleDeleteCollection: DeleteCollectionFn;
  isDeleting: boolean;
  isRevision: boolean;
  reorder: Reorder;
}

const MoreDropdown = ({
  collection,
  handleDeleteCollection,
  isDeleting,
  isRevision,
  reorder,
}: Props) => {
  const [anchorEl, setAnchorEl] = useState<null | HTMLButtonElement>(null);
  const open = Boolean(anchorEl);

  // Opens menu.
  const onOpen = (mouseEvent: MouseEvent<HTMLButtonElement>) => {
    setAnchorEl(mouseEvent.currentTarget);
  };

  // Closes menu.
  const onClose = () => {
    setAnchorEl(null);
  };

  return (
    <Fragment>
      <Button
        data-testid="collection-more-button"
        icon="DotsHorizontal"
        onClick={onOpen}
        open={open}
        sdsSize="small"
        sdsStyle="icon"
        sdsType="secondary"
      />
      <Menu
        collection={collection}
        handleDeleteCollection={handleDeleteCollection}
        isDeleting={isDeleting}
        isRevision={isRevision}
        menuProps={{ anchorEl, onClose, open }}
        reorder={reorder}
      />
    </Fragment>
  );
};

export default MoreDropdown;
