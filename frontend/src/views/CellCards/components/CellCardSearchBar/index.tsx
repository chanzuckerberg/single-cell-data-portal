import { Dropdown, DropdownPopper } from "czifui";
import { useState } from "react";

const OPTIONS = [
  {
    color: "#7057ff",
    description: "Good for newcomers",
    name: "good first issue",
  },
  {
    color: "#008672",
    description: "Extra attention is needed",
    name: "help wanted",
  },
];

export default function CellCardSearchBar(): JSX.Element {
  const [anchorEl, setAnchorEl] = useState<HTMLElement | null>(null);

  const openDropdown = (event: React.FocusEvent<HTMLElement>) => {
    console.log(event);
    setAnchorEl(event.currentTarget);
  };
  const closeDropdown = () => {};

  const handleSearch = (options: any) => {
    console.log(options);
  };
  const open = !!anchorEl;
  console.log(open);
  return (
    <div>
      <Dropdown
        DropdownMenuProps={{ loading: true, loadingText: "Loading ..." }}
        label="Click Target"
        onChange={() => {}}
        options={OPTIONS}
        buttonPosition="left"
        buttons
        closeOnBlur
        multiple
        search
        onClose={() => {}}
        InputDropdownProps={{ sdsStyle: "square" }}
        PopperComponent={DropdownPopper}
      />
    </div>
  );
}
