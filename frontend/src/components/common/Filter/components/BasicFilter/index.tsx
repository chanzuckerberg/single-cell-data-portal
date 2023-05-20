import { MouseEvent, ReactNode, useState } from "react";
import FilterLabel from "src/components/common/Filter/components/FilterLabel";
import { Filter, FilterPopover } from "../../common/style";

interface Props {
  content: ReactNode;
  isDisabled: boolean;
  label: string;
  tags: ReactNode;
  tooltip?: string;
}

export default function BasicFilter({
  content,
  isDisabled,
  label,
  tags,
  tooltip,
}: Props): JSX.Element {
  const [anchorEl, setAnchorEl] = useState<HTMLElement | null>(null);
  const isOpen = Boolean(anchorEl);
  return (
    <Filter>
      <FilterLabel
        isDisabled={isDisabled}
        isOpen={isOpen}
        label={label}
        onOpenFilter={(mouseEvent: MouseEvent<HTMLElement>) =>
          setAnchorEl(mouseEvent.currentTarget)
        }
        tooltip={tooltip}
      />
      <FilterPopover
        anchorEl={anchorEl}
        anchorOrigin={{
          horizontal: "right",
          vertical: "top",
        }}
        disableRestoreFocus
        disableScrollLock
        onClose={() => setAnchorEl(null)}
        open={isOpen}
        transformOrigin={{ horizontal: -5, vertical: 0 }}
      >
        {content}
      </FilterPopover>
      {tags}
    </Filter>
  );
}
