import React, { useMemo } from "react";

import { FilterOptionsState } from "@mui/material";
import { FilterOption } from "../../types";
import { CloseIcon, Tag, StyledTextField, StyledAutocomplete } from "./style";

interface Props {
  label: string;
  options: FilterOption[];
  selectedOptionIds: string[];
  handleChange: (options: FilterOption[]) => void;
}
function FilterDropdown({
  options,
  label,
  selectedOptionIds,
  handleChange,
}: Props): JSX.Element {
  const selectedOptions = useMemo(() => {
    return options.filter((option) => selectedOptionIds.includes(option.id));
  }, [options, selectedOptionIds]);

  return (
    <div>
      <StyledAutocomplete
        options={options}
        multiple
        onChange={(_: React.SyntheticEvent, newValue: FilterOption[]) => {
          handleChange(newValue);
        }}
        getOptionLabel={(option) => option.name}
        value={selectedOptions}
        disablePortal
        isOptionEqualToValue={(option, value) => option.id === value.id}
        renderInput={(params) => {
          return (
            <StyledTextField {...params} placeholder="Search" label={label} />
          );
        }}
        renderTags={(value, getTagProps) =>
          value.map((option, index) => (
            <Tag
              onClick={(event) => event.stopPropagation()}
              {...getTagProps({ index })}
              key={option.id}
            >
              {option.name}
              <CloseIcon
                onClick={() => {
                  handleChange(
                    selectedOptions.filter(
                      (selectedOption) => selectedOption.id !== option.id
                    )
                  );
                }}
              />
            </Tag>
          ))
        }
        filterOptions={(options, state) => {
          return options
            .filter((entity: FilterOption) => {
              const searchTerm = state.inputValue.toLowerCase();

              return (
                entity.name &&
                (entity.name.toLowerCase().includes(searchTerm) ||
                  entity.id.toLowerCase().includes(searchTerm))
              );
            })
            .sort((entityA: FilterOption, entityB: FilterOption) =>
              _sortOptions(entityA, entityB, state, selectedOptions)
            );
        }}
        disableCloseOnSelect
      />
    </div>
  );
}

function _sortOptions(
  entityA: FilterOption,
  entityB: FilterOption,
  state: FilterOptionsState<FilterOption>,
  selectedOptions: FilterOption[]
): number {
  const aRaw = entityA.name;
  const bRaw = entityB.name;
  const a = aRaw.toLowerCase();
  const b = bRaw.toLowerCase();
  const searchTerm = state.inputValue.toLowerCase();
  if (searchTerm === "") {
    // move selectedValues to top
    const includesA = selectedOptions
      .map((option) => option.name)
      .includes(aRaw);
    const includesB = selectedOptions
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
export default FilterDropdown;
