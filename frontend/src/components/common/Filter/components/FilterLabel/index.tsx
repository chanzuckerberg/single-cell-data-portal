import React, { MouseEvent } from "react";
import { ButtonDropdown } from "@czi-sds/components";
import FilterLabelTooltip from "src/components/common/Filter/components/FilterLabel/components/FilterLabelTooltip";
import { CategoryButton } from "./style";

type OnOpenFilterFn = (event: MouseEvent<HTMLElement>) => void;

interface Props {
  isDisabled: boolean;
  isOpen: boolean;
  label: string;
  onOpenFilter: OnOpenFilterFn;
  tooltip?: string;
}

export default function FilterLabel({
  isDisabled,
  isOpen,
  label,
  onOpenFilter,
  tooltip,
}: Props): JSX.Element {
  return (
    <FilterLabelTooltip tooltip={tooltip}>
      <CategoryButton isOpen={isOpen}>
        <ButtonDropdown
          disabled={isDisabled}
          onClick={onOpenFilter}
          sdsStyle="minimal"
          sdsType="secondary"
        >
          {label}
        </ButtonDropdown>
      </CategoryButton>
    </FilterLabelTooltip>
  );
}
