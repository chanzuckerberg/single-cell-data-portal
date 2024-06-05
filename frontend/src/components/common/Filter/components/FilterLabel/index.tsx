import React, { MouseEvent } from "react";
import FilterLabelTooltip from "src/components/common/Filter/components/FilterLabel/components/FilterLabelTooltip";
import { InputDropdown } from "src/components/common/Filter/components/FilterLabel/style";

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
      {/* The InputDropdown is enclosed within a <span> tag to enable tooltip functionality when the component is disabled. */}
      {/* See https://github.com/chanzuckerberg/sci-components/blob/main/packages/components/src/core/Tooltip/index.tsx#L28. */}
      <span>
        <InputDropdown
          disabled={isDisabled}
          intent="default"
          label={label}
          onClick={onOpenFilter}
          sdsStage="default"
          sdsStyle="minimal"
          sdsType="label"
          state="default"
        />
      </span>
    </FilterLabelTooltip>
  );
}
