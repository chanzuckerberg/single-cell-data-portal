import React, { useState } from "react";

import { Autocomplete, FilterOptionsState, TextField } from "@mui/material";
import { FilterOption } from "../../types";

const MAXIMUM_NUMBER_OF_SELECTED_OPTIONS = 2;
interface Props {
  label: string;
  options: FilterOption[];
  handleChange: (options: FilterOption[]) => void;
}
function FilterDropdown({ options, label, handleChange }: Props): JSX.Element {
  const [open, setOpen] = useState(false);
  const [selectedOptions, setSelectedOptions] = useState<FilterOption[]>([]);

  const handleFocus = () => {
    setOpen(true);
  };

  const handleBlur = () => {
    setOpen(false);
  };
  return (
    <div>
      <Autocomplete
        sx={{ width: 300 }}
        open={open}
        options={options}
        multiple
        onChange={(_: React.SyntheticEvent, newValue: FilterOption[]) => {
          setSelectedOptions(newValue);
          handleChange(newValue);
        }}
        getOptionLabel={(option) => option.name}
        value={selectedOptions}
        isOptionEqualToValue={(option, value) => option.id === value.id}
        renderInput={(params) => {
          const customEndAdornment = (
            <div
              onClick={() => {
                setOpen(!open);
              }}
            >
              {params.InputProps.endAdornment}
            </div>
          );
          return (
            <TextField
              {...params}
              onFocus={handleFocus}
              onBlur={handleBlur}
              InputProps={{
                ...params.InputProps,
                endAdornment: customEndAdornment,
              }}
              placeholder="Search"
              label={label}
            />
          );
        }}
        renderTags={(value, getTagProps) => [
          ...value
            .slice(0, MAXIMUM_NUMBER_OF_SELECTED_OPTIONS)
            .map((option, index) => (
              <div {...getTagProps({ index })} key={option.id}>
                {option.name}
              </div>
            )),
          ...(value.length > MAXIMUM_NUMBER_OF_SELECTED_OPTIONS
            ? [
                <div key="adornment">
                  +{value.length - MAXIMUM_NUMBER_OF_SELECTED_OPTIONS}
                </div>,
              ]
            : []),
        ]}
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
