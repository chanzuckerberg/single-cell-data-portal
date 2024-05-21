import { memo } from "react";

import { FlexRow, Wrapper } from "./style";
import { QueryGroup } from "src/views/DifferentialExpression/common/store/reducer";
import FilterDropdown from "./components/FilterDropdown";
import { Props } from "./types";
import { FilterOptionDimensions } from "../common/query_group_filter_dimensions";
import { QUERY_GROUP_KEY_TO_FILTER_DIMENSION_MAP } from "../common/constants";
import CopyButton from "./components/CopyButton";
import { DIFFERENTIAL_EXPRESSION_COPY_FILTERS_BUTTON_PREFIX } from "src/views/DifferentialExpression/common/constants";
import { QUERY_GROUP_KEYS, QUERY_GROUP_LABELS } from "./constants";
import { useConnect } from "./connect";

export default memo(function Filters({
  queryGroup,
  isQueryGroup1,
}: Props): JSX.Element {
  const { availableFilters, allAvailableFilters, handleFilterChange } =
    useConnect({ queryGroup, isQueryGroup1 });
  const components = QUERY_GROUP_KEYS.map((key, index) => {
    const queryGroupKey = key as keyof QueryGroup;
    const filterDimensionKey = QUERY_GROUP_KEY_TO_FILTER_DIMENSION_MAP[
      queryGroupKey
    ] as keyof FilterOptionDimensions;

    return {
      filterDropdownComponent: (
        <FilterDropdown
          key={`${queryGroupKey}-filter-dropdown`}
          label={QUERY_GROUP_LABELS[index]}
          options={availableFilters[filterDimensionKey]}
          allAvailableOptions={allAvailableFilters[filterDimensionKey]}
          selectedOptionIds={queryGroup[queryGroupKey]}
          handleChange={handleFilterChange(queryGroupKey)}
        />
      ),
      copyButtonComponent: isQueryGroup1 ? (
        <CopyButton
          key={`${queryGroupKey}-copy-button`}
          queryGroupKey={queryGroupKey}
          testId={`${DIFFERENTIAL_EXPRESSION_COPY_FILTERS_BUTTON_PREFIX}${queryGroupKey}`}
        />
      ) : null,
    };
  });

  return (
    <Wrapper>
      {components.map(({ filterDropdownComponent, copyButtonComponent }) => (
        <FlexRow key={filterDropdownComponent.key}>
          {isQueryGroup1 ? (
            <>
              {filterDropdownComponent}
              {copyButtonComponent}
            </>
          ) : (
            <>
              {copyButtonComponent}
              {filterDropdownComponent}
            </>
          )}
        </FlexRow>
      ))}
    </Wrapper>
  );
});
