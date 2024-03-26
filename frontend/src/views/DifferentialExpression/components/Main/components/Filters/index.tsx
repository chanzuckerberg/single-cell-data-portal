import { useState } from "react";
import isEqual from "lodash/isEqual";
import { memo, useCallback, useContext, useEffect, useMemo } from "react";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import {
  FilterDimensions,
  useQueryGroupFilterDimensions,
} from "src/common/queries/differentialExpression";
import { DispatchContext } from "src/views/DifferentialExpression/common/store";
import { Wrapper } from "./style";
import { QueryGroup } from "src/views/DifferentialExpression/common/store/reducer";
import {
  selectQueryGroup1Filters,
  selectQueryGroup2Filters,
} from "src/views/DifferentialExpression/common/store/actions";
import FilterDropdown from "./components/FilterDropdown";
import { FilterOption } from "./types";

const mapTermToFilterOption = (term: {
  id: string;
  name: string;
}): FilterOption => {
  return {
    name: term.name,
    id: term.id,
  };
};

type ModifiedFilterDimensions = Omit<
  FilterDimensions,
  "development_stage_terms" | "publication_citations"
> & {
  publication_citations: { name: string; id: string }[];
};

const EMPTY_OBJECT = {};

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

  const [availableFilters, setAvailableFilters] =
    useState<Partial<ModifiedFilterDimensions>>(EMPTY_OBJECT);

  const {
    data: {
      development_stage_terms: rawDevelopmentStages,
      disease_terms: rawDiseases,
      self_reported_ethnicity_terms: rawEthnicities,
      sex_terms: rawSexes,
      tissue_terms: rawTissues,
      cell_type_terms: rawCellTypes,
      publication_citations: rawPublicationCitations,
    },
    isLoading: rawIsLoading,
  } = useQueryGroupFilterDimensions(queryGroup);

  useEffect(() => {
    if (rawIsLoading) return;
    const newPublicationCitations = rawPublicationCitations.map((citation) => ({
      name: citation,
      id: citation,
    }));
    newPublicationCitations.sort((a, b) => a.name.localeCompare(b.name));

    const newSexes = rawSexes.map(mapTermToFilterOption);
    newSexes.sort((a, b) => a.name.localeCompare(b.name));

    const newDiseases = rawDiseases.map(mapTermToFilterOption);
    newDiseases.sort((a, b) =>
      a.name === "normal"
        ? -1
        : b.name === "normal"
        ? 1
        : a.name.localeCompare(b.name)
    );

    const newTissues = rawTissues.map(mapTermToFilterOption);
    newTissues.sort((a, b) => a.name.localeCompare(b.name));

    const newCellTypes = rawCellTypes.map(mapTermToFilterOption);
    newCellTypes.sort((a, b) => a.name.localeCompare(b.name));

    const newEthnicities = rawEthnicities.map(mapTermToFilterOption);
    newEthnicities.sort((a, b) => a.name.localeCompare(b.name));

    const newDevelopmentStages = rawDevelopmentStages;

    const newAvailableFilters = {
      development_stage_terms: newDevelopmentStages,
      disease_terms: newDiseases,
      self_reported_ethnicity_terms: newEthnicities,
      sex_terms: newSexes,
      tissue_terms: newTissues,
      cell_type_terms: newCellTypes,
      publication_citations: newPublicationCitations,
    };

    if (isEqual(availableFilters, newAvailableFilters)) return;

    setAvailableFilters(newAvailableFilters);
  }, [
    rawDiseases,
    rawEthnicities,
    rawSexes,
    rawCellTypes,
    rawIsLoading,
    rawDevelopmentStages,
    rawTissues,
    rawPublicationCitations,
    availableFilters,
    setAvailableFilters,
  ]);

  const {
    disease_terms = EMPTY_ARRAY,
    self_reported_ethnicity_terms = EMPTY_ARRAY,
    sex_terms = EMPTY_ARRAY,
    tissue_terms = EMPTY_ARRAY,
    cell_type_terms = EMPTY_ARRAY,
    publication_citations = EMPTY_ARRAY,
  } = availableFilters;

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
        const optionsWithNames = options.map((option) => {
          const typedOption = option as unknown as { id: string; name: string };
          const { id, name } = typedOption;
          return { id, name };
        });
        dispatch(selectQueryGroupFilters(key, optionsWithNames));
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
        handleChange={handleTissuesChange}
      />
      <FilterDropdown
        label="Cell Type"
        options={cell_type_terms}
        handleChange={handleCellTypesChange}
      />
      <FilterDropdown
        label="Publications"
        options={publication_citations}
        handleChange={handlePublicationCitationsChange}
      />
      <FilterDropdown
        label="Disease"
        options={disease_terms}
        handleChange={handleDiseasesChange}
      />
      <FilterDropdown
        label="Ethnicity"
        options={self_reported_ethnicity_terms}
        handleChange={handleEthnicitiesChange}
      />
      <FilterDropdown
        label="Sex"
        options={sex_terms}
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
