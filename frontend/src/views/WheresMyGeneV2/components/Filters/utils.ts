import { DefaultMenuSelectOption } from "@czi-sds/components";
import { FilterOption } from "./types";

export function isOptionEqualToValue(
  option: FilterOption,
  value: FilterOption
) {
  return option.id === value.id;
}
export function getOptionSelected(
  option: { id: string },
  value: { id: string }
): boolean {
  return option.id === value.id;
}

export const mapTermToFilterOption = (term: {
  id: string;
  name: string;
}): FilterOption => {
  return {
    name: term.name,
    label: `${term.name} (${term.id})`,
    id: term.id,
  };
};

export function sortOptions(
  a: DefaultMenuSelectOption,
  b: DefaultMenuSelectOption
) {
  if (a.name < b.name) {
    return -1;
  }
  if (a.name > b.name) {
    return 1;
  }
  return 0;
}
