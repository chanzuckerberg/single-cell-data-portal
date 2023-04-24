import { DefaultDropdownMenuOption, InputSearch } from "czifui";
import { useState, useEffect, useRef } from "react";
import { noop } from "src/common/constants/utils";
import { SectionItem, SectionTitle, StyledPopper } from "./style";
import { allCellTypes } from "./fixture";

export default function CellCardSearchBar(): JSX.Element {
  const dropdownRef = useRef(null);
  const [anchorEl, setAnchorEl] = useState<HTMLElement | null>(
    dropdownRef.current
  );
  const [open, setOpen] = useState(false);

  const [value, setValue] = useState("");

  const handleFocus = (event: React.FocusEvent<HTMLInputElement>) => {
    if (!anchorEl)
      setAnchorEl(event.currentTarget.parentElement?.parentElement ?? null);
    setOpen(true);
  };
  const handleBlur = () => {
    setTimeout(() => {
      setOpen(false);
    }, 100);
  };
  const handleEscape = (event: React.KeyboardEvent<HTMLInputElement>) => {
    if (event.key === "Escape") {
      setOpen(false);
      (event.target as HTMLElement).blur();
    }
  };
  const handleSearch = (event: React.ChangeEvent<HTMLInputElement>) => {
    setValue(event.target.value);
  };

  const handleClick = (event: React.MouseEvent<HTMLDivElement, MouseEvent>) => {
    console.log((event.target as HTMLDivElement).innerText);
  };

  return (
    <div>
      <InputSearch
        sdsStyle="square"
        placeholder="Search cell type or tissues"
        variant="outlined"
        id="cell-cards-search-bar"
        fullWidth
        label="Cell Cards"
        onChange={handleSearch}
        onFocus={handleFocus}
        onBlur={handleBlur}
        onKeyDown={handleEscape}
        autoComplete="off"
      />
      <StyledPopper
        disablePortal={true}
        ref={dropdownRef}
        open={open}
        anchorEl={anchorEl}
        onResize={noop}
        onResizeCapture={noop}
        width={anchorEl?.clientWidth || 0}
        placement="bottom-start"
        modifiers={[
          {
            name: "flip",
            enabled: false,
          },
        ]}
      >
        <SectionTitle>Cell Types</SectionTitle>
        {allCellTypes
          .filter((cellType) => {
            const [cellTypeName] = Object.values(cellType);
            return cellTypeName.toLowerCase().includes(value.toLowerCase());
          })
          .map((cellType) => {
            const [cellTypeId] = Object.keys(cellType);
            const [cellTypeName] = Object.values(cellType);
            return (
              <SectionItem id={cellTypeId} onClick={handleClick}>
                {cellTypeName}
              </SectionItem>
            );
          })}
      </StyledPopper>
    </div>
  );
}
