import { Popover, Position } from "@blueprintjs/core";
import { ReactNode } from "react";
import { Filter } from "../../common/style";

interface Props {
  content: ReactNode;
  isDisabled: boolean;
  tags: ReactNode;
  target: ReactNode;
}

export default function BasicFilter({
  content,
  isDisabled,
  tags,
  target,
}: Props): JSX.Element {
  return (
    <Filter>
      <Popover
        boundary="viewport"
        disabled={isDisabled}
        minimal
        modifiers={{ offset: { offset: "0, 5" } }}
        position={Position.RIGHT}
      >
        {target}
        {content}
      </Popover>
      {tags}
    </Filter>
  );
}
