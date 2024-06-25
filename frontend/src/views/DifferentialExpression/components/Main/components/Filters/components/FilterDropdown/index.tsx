import React from "react";

import { FilterOption } from "../../types";
import {
  CloseIcon,
  PrimaryTag,
  GrayTag,
  StyledTextField,
  StyledAutocomplete,
  StyledCallout,
} from "./style";
import {
  DIFFERENTIAL_EXPRESSION_FILTER_AUTOCOMPLETE_PREFIX,
  DIFFERENTIAL_EXPRESSION_FILTER_TAG_GRAY,
  DIFFERENTIAL_EXPRESSION_FILTER_TAG_PRIMARY,
} from "src/views/DifferentialExpression/common/constants";
import { sortOptions } from "./utils";
import { useConnect } from "./connect";
import { Props } from "./types";
import { Icon } from "@czi-sds/components";

function FilterDropdown({
  options,
  label,
  allAvailableOptions,
  selectedOptionIds,
  handleChange,
}: Props): JSX.Element {
  const {
    selectedOptions,
    previousSelectedOptions,
    setPreviousSelectedOptions,
  } = useConnect({
    options,
    allAvailableOptions,
    selectedOptionIds,
  });

  return (
    <div>
      <StyledAutocomplete
        data-testid={`${DIFFERENTIAL_EXPRESSION_FILTER_AUTOCOMPLETE_PREFIX}${label}`}
        options={options}
        onClose={() => setPreviousSelectedOptions(selectedOptions)}
        multiple
        onChange={(_: React.SyntheticEvent, newValue: FilterOption[]) => {
          handleChange(newValue);
          if (newValue.length === 0) {
            setPreviousSelectedOptions(newValue);
          }
        }}
        getOptionLabel={(option) => option.name}
        value={selectedOptions}
        disablePortal
        noOptionsText={
          <StyledCallout
            icon={
              <Icon
                color="gray"
                sdsIcon="InfoCircle"
                sdsSize="l"
                sdsType="static"
              />
            }
            intent="info"
          >
            Results may be limited by other selections made in this cell group.
          </StyledCallout>
        }
        isOptionEqualToValue={(option, value) => option.id === value.id}
        popupIcon={null}
        renderInput={(params) => {
          return (
            <StyledTextField {...params} placeholder="Search" label={label} />
          );
        }}
        renderTags={(value, getTagProps) =>
          value.map((option, index) => {
            const isUnavailable = option.unavailable;
            const TagComponent = isUnavailable ? GrayTag : PrimaryTag;
            return (
              <TagComponent
                data-testid={
                  isUnavailable
                    ? DIFFERENTIAL_EXPRESSION_FILTER_TAG_GRAY
                    : DIFFERENTIAL_EXPRESSION_FILTER_TAG_PRIMARY
                }
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
              </TagComponent>
            );
          })
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
              sortOptions(entityA, entityB, state, previousSelectedOptions)
            );
        }}
        disableCloseOnSelect
      />
    </div>
  );
}

export default FilterDropdown;
