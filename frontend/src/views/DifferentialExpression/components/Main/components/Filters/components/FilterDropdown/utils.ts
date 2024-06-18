import { FilterOption } from "../../types";
import { FilterOptionsState } from "@mui/material";

export function sortOptions(
  entityA: FilterOption,
  entityB: FilterOption,
  state: FilterOptionsState<FilterOption>,
  previousSelectedOptions: FilterOption[]
): number {
  const aRaw = entityA.name;
  const bRaw = entityB.name;
  const a = aRaw.toLowerCase();
  const b = bRaw.toLowerCase();
  const searchTerm = state.inputValue.toLowerCase();
  if (searchTerm === "") {
    // move selectedValues to top
    const includesA = previousSelectedOptions
      .map((option) => option.name)
      .includes(aRaw);
    const includesB = previousSelectedOptions
      .map((option) => option.name)
      .includes(bRaw);
    if (includesA && includesB) {
      return a.localeCompare(b);
    }
    if (includesA) {
      return -1;
    }
    if (includesB) {
      return 1;
    }
  }
  // Determine if each item starts with the search term
  const aStartsWithSearch = a.startsWith(searchTerm);
  const bStartsWithSearch = b.startsWith(searchTerm);

  if (aStartsWithSearch && !bStartsWithSearch) {
    return -1;
  }
  if (!aStartsWithSearch && bStartsWithSearch) {
    return 1;
  }
  return a.localeCompare(b);
}
