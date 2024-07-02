import { Fragment } from "react";
import { Props } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/MoreDropdown/types";
import {
  DATASET_MORE_BUTTON,
  MORE_ICON_BUTTON_PROPS,
} from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/MoreDropdown/constants";
import { StyledButton } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/MoreDropdown/style";
import { useMenu } from "src/views/Collection/hooks/useMenu";
import Menu from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/MoreDropdown/components/Menu";

export default function MoreDropdown({
  collectionId,
  dataset,
  disabled,
  menuItemProps,
}: Props): JSX.Element {
  const { anchorEl, onClose, onOpen, open } = useMenu();
  return (
    <Fragment>
      <StyledButton
        {...MORE_ICON_BUTTON_PROPS}
        data-testid={DATASET_MORE_BUTTON}
        disabled={disabled}
        onClick={onOpen}
        open={open}
      />
      <Menu
        collectionId={collectionId}
        dataset={dataset}
        menuItemProps={menuItemProps}
        menuProps={{ anchorEl, onClose, open }}
      />
    </Fragment>
  );
}
