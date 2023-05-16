import { useRef, useState } from "react";
import { TextField } from "@mui/material";
import SearchIcon from "@mui/icons-material/Search";
import InputAdornment from "@mui/material/InputAdornment";
import { SectionItem, StyledAutocomplete } from "./style";
import { ROUTES } from "src/common/constants/routes";
import { useCellTypes } from "src/common/queries/cellCards";
import { useRouter } from "next/router";

interface CellType {
  id: string;
  label: string;
}

export default function CellCardSearchBar(): JSX.Element {
  const router = useRouter();
  const { data: cellTypes } = useCellTypes();

  const [open, setOpen] = useState(false);

  // Used for keyboard navigation for cell type search
  const [highlightedCellType, setHighlightedCellType] =
    useState<CellType | null>(null);

  const handleFocus = () => {
    setOpen(true);
  };

  const handleBlur = () => {
    setOpen(false);
  };

  function changeCellType(cellTypeId: string) {
    if (cellTypeId) {
      router.push(`${ROUTES.CELL_CARDS}/${cellTypeId.replace(":", "_")}`);
      document.getElementById("cell-cards-search-bar")?.blur();
    }
  }

  return (
    <div>
      <StyledAutocomplete
        open={open}
        onKeyDown={(event) => {
          if (highlightedCellType && event.key === "Enter") {
            changeCellType(highlightedCellType.id);
          }
        }}
        onHighlightChange={(_, value) => {
          const cellType = value as CellType;
          setHighlightedCellType(cellType);
        }}
        disablePortal
        id="cell-cards-search-bar"
        options={cellTypes ?? []}
        renderInput={(params) => (
          <TextField
            {...params}
            onFocus={handleFocus}
            onBlur={handleBlur}
            InputProps={{
              ...params.InputProps,
              endAdornment: (
                <InputAdornment position="end">
                  <SearchIcon />
                </InputAdornment>
              ),
            }}
            label="Search cell types or tissues"
          />
        )}
        renderOption={(props, option) => {
          const cellType = option as CellType;
          return (
            <SectionItem
              {...props}
              key={cellType.id}
              onClick={() => {
                changeCellType(cellType.id);
              }}
            >
              {cellType.label}
            </SectionItem>
          );
        }}
        autoComplete
        filterOptions={(options, state) => {
          return options
            .filter((option) => {
              const cellType = option as CellType;
              return (
                cellType.label &&
                cellType.label
                  .toLowerCase()
                  .includes(state.inputValue.toLowerCase())
              );
            })
            .sort((cellTypeA, cellTypeB) => {
              const aRaw = (cellTypeA as CellType).label;
              const bRaw = (cellTypeB as CellType).label;
              const a = aRaw.toLowerCase();
              const b = bRaw.toLowerCase();
              const searchTerm = state.inputValue.toLowerCase();
              if (a.startsWith(searchTerm) && !b.startsWith(searchTerm)) {
                return -1;
              }
              if (!a.startsWith(searchTerm) && b.startsWith(searchTerm)) {
                return 1;
              }
              return a.localeCompare(b);
            });
        }}
      />
    </div>
  );
}
