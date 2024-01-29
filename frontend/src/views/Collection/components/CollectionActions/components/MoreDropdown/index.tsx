import { Fragment, MouseEvent, useState } from "react";
import { Collection } from "src/common/entities";
import Menu from "./components/Menu";
import { DeleteCollectionFn } from "src/views/Collection/components/CollectionActions";
import { ButtonIcon } from "src/views/Collection/components/CollectionActions/components/MoreDropdown/style";
import { OnSetReorderModeFn } from "src/common/hooks/useReorderMode";

interface Props {
  collection: Collection;
  handleDeleteCollection: DeleteCollectionFn;
  isDeleting: boolean;
  isReorderUX: boolean;
  isRevision: boolean;
  onSetReorderMode: OnSetReorderModeFn;
}

const MoreDropdown = ({
  collection,
  handleDeleteCollection,
  isDeleting,
  isReorderUX,
  isRevision,
  onSetReorderMode,
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
      <ButtonIcon
        data-testid="collection-more-button"
        onClick={onOpen}
        open={open}
        sdsIcon="dotsHorizontal"
        sdsSize="small"
        sdsType="tertiary"
      />
      <Menu
        collection={collection}
        handleDeleteCollection={handleDeleteCollection}
        isDeleting={isDeleting}
        isReorderUX={isReorderUX}
        isRevision={isRevision}
        menuProps={{ anchorEl, onClose, open }}
        onSetReorderMode={onSetReorderMode}
      />
    </Fragment>
  );
};

export default MoreDropdown;
