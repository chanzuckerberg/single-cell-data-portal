import { Popover, Position } from "@blueprintjs/core";
import { ElementType, ReactNode } from "react";
import { useFilterSearch } from "src/components/common/Filter/components/FilterSearch/common/useFilterSearch";
import { Filter } from "../../common/style";

interface Props {
  Content: ElementType;
  flipEnabled?: boolean;
  isDisabled: boolean;
  tags: ReactNode;
  target: ReactNode;
}

export default function BasicFilter({
  Content,
  flipEnabled = true,
  isDisabled,
  tags,
  target,
}: Props): JSX.Element {
  const filterSearchState = useFilterSearch();
  const { clearSearchValueFn } = filterSearchState;
  return (
    <Filter>
      <Popover
        boundary="viewport"
        disabled={isDisabled}
        minimal
        modifiers={{
          flip: { enabled: flipEnabled },
          offset: { offset: "0, 4" },
        }}
        onClosed={clearSearchValueFn}
        position={Position.BOTTOM_LEFT}
      >
        {target}
        <Content {...filterSearchState} />
      </Popover>
      {tags}
    </Filter>
  );
}
