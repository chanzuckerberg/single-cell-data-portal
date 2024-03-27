import { memo, useCallback, useContext, useMemo } from "react";

import { DispatchContext } from "src/views/DifferentialExpression/common/store";
import { Wrapper } from "./style";
import { QueryGroup } from "src/views/DifferentialExpression/common/store/reducer";
import {
  selectQueryGroup1Filters,
  selectQueryGroup2Filters,
} from "src/views/DifferentialExpression/common/store/actions";
import FilterDropdown from "./components/FilterDropdown";
import { FilterOption } from "./types";
import useProcessedQueryGroupFilterDimensions from "../common/query_group_filter_dimensions";

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

  const {
    disease_terms,
    self_reported_ethnicity_terms,
    sex_terms,
    tissue_terms,
    cell_type_terms,
    publication_citations,
  } = useProcessedQueryGroupFilterDimensions(queryGroup);

  const handleFilterChange = useCallback(
    function handleFilterChange_(
      key: keyof QueryGroup
    ): (options: FilterOption[] | undefined) => void {
      let currentOptions: FilterOption[] | undefined;

      return (options: FilterOption[] | undefined): void => {
        if (
          !dispatch ||
          !options ||
          // If the options are the same
          JSON.stringify(options.sort(sortOptions)) ===
            JSON.stringify(currentOptions?.sort(sortOptions)) ||
          // If the options change from null to [], which is the default value
          (!currentOptions && JSON.stringify(options) === "[]")
        ) {
          return;
        }

        currentOptions = options;

        dispatch(selectQueryGroupFilters(key, currentOptions));
      };
    },
    [dispatch, selectQueryGroupFilters]
  );

  const handlePublicationCitationsChange = useMemo(
    () => handleFilterChange("publicationCitations"),
    [handleFilterChange]
  );

  const handleDiseasesChange = useMemo(
    () => handleFilterChange("diseases"),
    [handleFilterChange]
  );

  const handleEthnicitiesChange = useMemo(
    () => handleFilterChange("ethnicities"),
    [handleFilterChange]
  );

  const handleSexesChange = useMemo(
    () => handleFilterChange("sexes"),
    [handleFilterChange]
  );

  const handleTissuesChange = useMemo(
    () => handleFilterChange("tissues"),
    [handleFilterChange]
  );

  const handleCellTypesChange = useMemo(
    () => handleFilterChange("cellTypes"),
    [handleFilterChange]
  );

  return (
    <Wrapper>
      <FilterDropdown
        label="Tissue"
        options={tissue_terms}
        selectedOptionIds={queryGroup.tissues}
        handleChange={handleTissuesChange}
      />
      <FilterDropdown
        label="Cell Type"
        options={cell_type_terms}
        selectedOptionIds={queryGroup.cellTypes}
        handleChange={handleCellTypesChange}
      />
      <FilterDropdown
        label="Publications"
        options={publication_citations}
        selectedOptionIds={queryGroup.publicationCitations}
        handleChange={handlePublicationCitationsChange}
      />
      <FilterDropdown
        label="Disease"
        options={disease_terms}
        selectedOptionIds={queryGroup.diseases}
        handleChange={handleDiseasesChange}
      />
      <FilterDropdown
        label="Ethnicity"
        options={self_reported_ethnicity_terms}
        selectedOptionIds={queryGroup.ethnicities}
        handleChange={handleEthnicitiesChange}
      />
      <FilterDropdown
        label="Sex"
        options={sex_terms}
        selectedOptionIds={queryGroup.sexes}
        handleChange={handleSexesChange}
      />
    </Wrapper>
  );
});

function sortOptions(a: FilterOption, b: FilterOption) {
  if (a.name < b.name) {
    return -1;
  }
  if (a.name > b.name) {
    return 1;
  }
  return 0;
}
