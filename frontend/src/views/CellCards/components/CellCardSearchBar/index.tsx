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

  return (
    <div>
      <StyledAutocomplete
        disablePortal
        id="cell-cards-search-bar"
        options={cellTypes ?? []}
        renderInput={(params) => (
          <TextField
            {...params}
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
        renderOption={(_, option) => {
          const cellType = option as CellType;
          return (
            <SectionItem
              onClick={() => {
                router.push(
                  `${ROUTES.CELL_CARDS}/${cellType.id.replace(":", "_")}`
                );
                document.getElementById("cell-cards-search-bar")?.blur();
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
              return (option as CellType).label
                .toLowerCase()
                .includes(state.inputValue.toLowerCase());
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
