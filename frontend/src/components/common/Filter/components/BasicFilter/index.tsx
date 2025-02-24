import {
  MouseEvent,
  ReactNode,
  useEffect,
  useState,
  useRef,
  useContext,
} from "react";
import FilterLabel from "src/components/common/Filter/components/FilterLabel";
import { Filter, FilterPopover } from "../../common/style";
import { FilterWrapper } from "./style";
import { FilterControlContext } from "src/components/common/Filter/context";

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
  const { openSpecificFilter, setOpenSpecificFilter } =
    useContext(FilterControlContext);

  const anchorElRef = useRef<HTMLDivElement | null>(null);

  useEffect(() => {
    if (openSpecificFilter && label !== openSpecificFilter) {
      setAnchorEl(null);
    } else if (openSpecificFilter === label) {
      setAnchorEl(anchorElRef.current);
    }
  }, [openSpecificFilter, label]);

  const onClose = () => {
    setOpenSpecificFilter(null);
    setAnchorEl(null);
  };

  const isOpen = Boolean(anchorEl);
  return (
    <Filter>
      <FilterWrapper ref={anchorElRef}>
        <FilterLabel
          isDisabled={isDisabled}
          label={label}
          onOpenFilter={(mouseEvent: MouseEvent<HTMLElement>) =>
            setAnchorEl(mouseEvent.currentTarget)
          }
          tooltip={tooltip}
        />
      </FilterWrapper>
      <FilterPopover
        anchorEl={anchorEl}
        anchorOrigin={{
          horizontal: "right",
          vertical: "top",
        }}
        disableRestoreFocus
        disableScrollLock
        onClose={onClose}
        open={isOpen}
        transformOrigin={{ horizontal: -5, vertical: 0 }}
      >
        {content}
      </FilterPopover>
      {tags}
    </Filter>
  );
}
