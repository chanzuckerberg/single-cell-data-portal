import { memo, useCallback, useContext } from "react";

import { DispatchContext } from "src/views/DifferentialExpression/common/store";
import { FlexRow, Wrapper } from "./style";
import {
  EMPTY_FILTERS,
  QueryGroup,
} from "src/views/DifferentialExpression/common/store/reducer";
import {
  selectQueryGroup1Filters,
  selectQueryGroup2Filters,
} from "src/views/DifferentialExpression/common/store/actions";
import FilterDropdown from "./components/FilterDropdown";
import { FilterOption } from "./types";
import useProcessedQueryGroupFilterDimensions, {
  FilterOptionDimensions,
} from "../common/query_group_filter_dimensions";
import { QUERY_GROUP_KEY_TO_FILTER_DIMENSION_MAP } from "../common/constants";
import CopyInvertButton from "./components/CopyInvertButton";

const QUERY_GROUP_KEYS = [
  "tissues",
  "cellTypes",
  "publicationCitations",
  "diseases",
  "ethnicities",
  "sexes",
];
const QUERY_GROUP_LABELS = [
  "Tissue",
  "Cell Type",
  "Publications",
  "Disease",
  "Ethnicity",
  "Sex",
];

interface Props {
  queryGroup: QueryGroup;
  isQueryGroup1: boolean;
}
export default memo(function Filters({
  queryGroup,
  isQueryGroup1,
}: Props): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const selectQueryGroupFilters = isQueryGroup1
    ? selectQueryGroup1Filters
    : selectQueryGroup2Filters;

  const { availableFilters } =
    useProcessedQueryGroupFilterDimensions(queryGroup);

  const handleFilterChange = useCallback(
    function handleFilterChange_(
      key: keyof QueryGroup
    ): (options: FilterOption[] | undefined) => void {
      return (options: FilterOption[] | undefined): void => {
        if (!dispatch || !options) {
          return;
        }

        dispatch(selectQueryGroupFilters(key, options));
      };
    },
    [dispatch, selectQueryGroupFilters]
  );

  const { availableFilters: allAvailableFilters } =
    useProcessedQueryGroupFilterDimensions(EMPTY_FILTERS);

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
      copyInvertButtonComponent: (
        <CopyInvertButton
          key={`${queryGroupKey}-copy-invert-button`}
          queryGroupKey={queryGroupKey}
          isCopy={isQueryGroup1}
        />
      ),
    };
  });

  return (
    <Wrapper>
      {components.map(
        ({ filterDropdownComponent, copyInvertButtonComponent }) => (
          <FlexRow key={filterDropdownComponent.key}>
            {isQueryGroup1 ? (
              <>
                {filterDropdownComponent}
                {copyInvertButtonComponent}
              </>
            ) : (
              <>
                {copyInvertButtonComponent}
                {filterDropdownComponent}
              </>
            )}
          </FlexRow>
        )
      )}
    </Wrapper>
  );
});
