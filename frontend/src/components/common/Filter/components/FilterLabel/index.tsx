import { Button } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { CATEGORY_LABEL } from "src/components/common/Filter/common/entities";
import FilterLabelTooltip from "src/components/common/Filter/components/FilterLabel/FilterLabelTooltip";
import { CategoryButton } from "./style";

interface Props {
  isDisabled: boolean;
  label: CATEGORY_LABEL;
  tooltip?: string;
}

export default function FilterLabel({
  isDisabled,
  label,
  tooltip,
}: Props): JSX.Element {
  return (
    <FilterLabelTooltip tooltip={tooltip}>
      <CategoryButton>
        <Button
          disabled={isDisabled}
          minimal
          rightIcon={IconNames.CHEVRON_DOWN}
          text={label}
        />
      </CategoryButton>
    </FilterLabelTooltip>
  );
}
