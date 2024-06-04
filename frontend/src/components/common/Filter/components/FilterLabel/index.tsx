import React, { MouseEvent } from "react";
import FilterLabelTooltip from "src/components/common/Filter/components/FilterLabel/components/FilterLabelTooltip";
import { ButtonDropdown } from "./style";

type OnOpenFilterFn = (event: MouseEvent<HTMLElement>) => void;

interface Props {
  isDisabled: boolean;
  label: string;
  onOpenFilter: OnOpenFilterFn;
  tooltip?: string;
}

export default function FilterLabel({
  isDisabled,
  label,
  onOpenFilter,
  tooltip,
}: Props): JSX.Element {
  return (
    <FilterLabelTooltip tooltip={tooltip}>
      <ButtonDropdown
        disabled={isDisabled}
        onClick={onOpenFilter}
        sdsStyle="minimal"
        sdsType="secondary"
      >
        {label}
      </ButtonDropdown>
    </FilterLabelTooltip>
  );
}
