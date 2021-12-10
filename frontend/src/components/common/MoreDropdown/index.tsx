import {
  Button,
  IButtonProps,
  IPopoverProps,
  Popover,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { FEATURES } from "src/common/featureFlags/features";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { ActionButton } from "src/components/common/Grid/components/ActionButton/style";

interface Props {
  popoverProps?: IPopoverProps;
  buttonProps?: IButtonProps;
}

const MoreDropdown = ({
  popoverProps = {},
  buttonProps = {},
}: Props): JSX.Element => {
  const isFilterEnabled = useFeatureFlag(FEATURES.FILTER);
  const MoreButton = isFilterEnabled ? ActionButton : Button;
  return (
    <Popover {...popoverProps}>
      <MoreButton
        {...buttonProps}
        minimal
        icon={IconNames.MORE}
        data-test-id="collection-more-button"
      />
    </Popover>
  );
};

export default MoreDropdown;
