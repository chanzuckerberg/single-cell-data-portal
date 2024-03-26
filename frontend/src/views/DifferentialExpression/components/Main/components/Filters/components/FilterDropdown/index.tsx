import React, { useState } from "react";

import { Autocomplete, TextField } from "@mui/material";
import { FilterOption } from "../../types";

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
        sx={{ width: 300, height: 300 }}
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
          const customLabel = label;
          return (
            <TextField
              {...params}
              sx={{ width: 300, height: 300 }}
              onFocus={handleFocus}
              onBlur={handleBlur}
              InputProps={{
                ...params.InputProps,
                endAdornment: customEndAdornment,
              }}
              placeholder="Search"
              label={customLabel}
            />
          );
        }}
        renderTags={(value, getTagProps) =>
          value.map((option, index) => (
            <div {...getTagProps({ index })} key={option.id}>
              {option.name}
            </div>
          ))
        }
      />
    </div>
  );
}

export default FilterDropdown;
