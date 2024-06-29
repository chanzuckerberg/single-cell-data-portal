import { Fragment } from "react";
import { Props } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/MoreDropdown/types";
import {
  DATASET_MORE_BUTTON,
  DEFAULT_BUTTON_PROPS,
} from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/MoreDropdown/constants";
import { StyledButton } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/MoreDropdown/style";
import { useMoreMenu } from "src/views/Collection/hooks/useMoreMenu";
import Menu from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/MoreDropdown/components/Menu";

export default function MoreDropdown({
  collectionId,
  dataset,
  disabled,
  menuItemProps,
}: Props): JSX.Element {
  const { anchorEl, onClose, onOpen, open } = useMoreMenu();
  return (
    <Fragment>
      <StyledButton
        {...DEFAULT_BUTTON_PROPS}
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
