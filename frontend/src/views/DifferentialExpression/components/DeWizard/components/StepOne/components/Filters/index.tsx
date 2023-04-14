import { useState } from "react";
import { createFilterOptions } from "@mui/material";
import {
  ComplexFilterInputDropdown,
  DefaultMenuSelectOption,
  InputDropdownProps,
} from "czifui";
import isEqual from "lodash/isEqual";
import { memo, useCallback, useContext, useEffect, useMemo } from "react";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import {
  FilterDimensions,
  RawDataset,
  useFilterDimensions,
} from "src/common/queries/differentialExpression";
import {
  DispatchContext,
  StateContext,
} from "src/views/DifferentialExpression/common/store";
import { selectFilters } from "src/views/DifferentialExpression/common/store/actions";
import { Filters as IFilters } from "src/views/DifferentialExpression/common/types";
import Organism from "./components/Organism";
import {
  StyledComplexFilter,
  StyledComplexFilterInputDropdown,
  Wrapper,
  StyledPopper,
  FilterLabel,
} from "./style";

const filterOptions = createFilterOptions({
  stringify: (option: RawDataset) =>
    `${option.label} ${option.collection_label}`,
});

const DropdownMenuProps = {
  filterOptions,
  getOptionSelected,
};

interface FilterOption {
  name: string;
  label: string;
  id: string;
}

const mapTermToFilterOption = (term: {
  id: string;
  name: string;
}): FilterOption => {
  return {
    name: term.name,
    label: `${term.name} (${term.id})`,
    id: term.id,
  };
};

const EMPTY_OBJECT = {};

export default memo(function Filters(): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const state = useContext(StateContext);

  const [availableFilters, setAvailableFilters] =
    useState<Partial<FilterDimensions>>(EMPTY_OBJECT);

  const { selectedFilters } = state;

  const {
    datasets: datasetIds,
    diseases,
    ethnicities,
    sexes,
    tissues,
    developmentStages,
  } = selectedFilters;

  const {
    data: {
      datasets: rawDatasets,
      development_stage_terms: rawDevelopmentStages,
      disease_terms: rawDiseases,
      self_reported_ethnicity_terms: rawEthnicities,
      sex_terms: rawSexes,
      tissue_terms: rawTissues,
    },
    isLoading: rawIsLoading,
  } = useFilterDimensions();

  const InputDropdownProps = {
    sdsStyle: "minimal",
  } as Partial<InputDropdownProps>;

  // (thuang): We only update available filters when API call is done,
  // otherwise when `useFilterDimensions()` is still loading, its filters
  // will temporarily be empty, and thus resetting the selected filter values
  useEffect(() => {
    if (rawIsLoading) return;
    const newDatasets = rawDatasets.map((dataset) => ({
      ...dataset,
      details: dataset.collection_label,
      name: dataset.label,
    }));
    newDatasets.sort((a, b) => a.name.localeCompare(b.name));

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

    const newEthnicities = rawEthnicities.map(mapTermToFilterOption);
    newEthnicities.sort((a, b) => a.name.localeCompare(b.name));

    const newDevelopmentStages = rawDevelopmentStages.map(
      mapTermToFilterOption
    );
    newDevelopmentStages.sort((a, b) => a.name.localeCompare(b.name));

    const newAvailableFilters = {
      datasets: newDatasets,
      development_stage_terms: newDevelopmentStages,
      disease_terms: newDiseases,
      self_reported_ethnicity_terms: newEthnicities,
      sex_terms: newSexes,
      tissue_terms: newTissues,
    };

    if (isEqual(availableFilters, newAvailableFilters)) return;

    setAvailableFilters(newAvailableFilters);
  }, [rawDatasets, rawDevelopmentStages, rawDiseases, rawEthnicities, rawSexes, rawIsLoading, availableFilters, setAvailableFilters]);

  const {
    datasets = EMPTY_ARRAY,
    disease_terms = EMPTY_ARRAY,
    self_reported_ethnicity_terms = EMPTY_ARRAY,
    sex_terms = EMPTY_ARRAY,
    tissue_terms = EMPTY_ARRAY,
    development_stage_terms = EMPTY_ARRAY,
  } = availableFilters;

  const selectedDatasets = useMemo(() => {
    return datasets.filter((dataset) => datasetIds?.includes(dataset.id));
  }, [datasets, datasetIds]);

  const selectedDiseases = useMemo(() => {
    return disease_terms.filter((disease) => diseases?.includes(disease.id));
  }, [disease_terms, diseases]);

  const selectedEthnicities = useMemo(() => {
    return self_reported_ethnicity_terms.filter((ethnicity) =>
      ethnicities?.includes(ethnicity.id)
    );
  }, [self_reported_ethnicity_terms, ethnicities]);

  const selectedSexes = useMemo(() => {
    return sex_terms.filter((sex) => sexes?.includes(sex.id));
  }, [sex_terms, sexes]);

  const selectedDevelopmentStages = useMemo(() => {
    return development_stage_terms.filter((stage) =>
      developmentStages?.includes(stage.id)
    );
  }, [development_stage_terms, developmentStages]);

  const selectedTissues = useMemo(() => {
    return tissue_terms.filter((tissue) => tissues?.includes(tissue.id));
  }, [tissue_terms, tissues]);

  const handleFilterChange = useCallback(
    function handleFilterChange_(
      key: keyof IFilters
    ): (options: DefaultMenuSelectOption[] | null) => void {
      let currentOptions: DefaultMenuSelectOption[] | null = null;

      return (options: DefaultMenuSelectOption[] | null): void => {
        if (
          !dispatch ||
          !options ||
          // If the options are the same
          JSON.stringify(options.sort(sortOptions)) ===
            JSON.stringify(currentOptions?.sort(sortOptions)) ||
          // If the options change from null to [], which is the default value
          (currentOptions === null && JSON.stringify(options) === "[]")
        ) {
          return;
        }

        currentOptions = options;

        dispatch(
          selectFilters(
            key,
            options.map((option) => (option as unknown as { id: string }).id)
          )
        );
      };
    },
    [dispatch]
  );

  const handleDatasetsChange = useMemo(
    () => handleFilterChange("datasets"),
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

  const handleDevelopmentStagesChange = useMemo(
    () => handleFilterChange("developmentStages"),
    [handleFilterChange]
  );

  return (
    <Wrapper>
      <Organism />
      <div>
        <FilterLabel> Filters </FilterLabel>
        <StyledComplexFilter
          multiple
          data-testid="de-tissue-filter"
          search
          label="Tissue"
          options={tissue_terms as unknown as DefaultMenuSelectOption[]}
          onChange={handleTissuesChange}
          value={selectedTissues as unknown as DefaultMenuSelectOption[]}
          InputDropdownComponent={
            StyledComplexFilterInputDropdown as typeof ComplexFilterInputDropdown
          }
          DropdownMenuProps={DropdownMenuProps}
          InputDropdownProps={InputDropdownProps}
          PopperComponent={StyledPopper}
        />
        <StyledComplexFilter
          multiple
          data-testid="de-dataset-filter"
          search
          label="Dataset"
          options={datasets as unknown as DefaultMenuSelectOption[]}
          onChange={handleDatasetsChange}
          value={selectedDatasets as unknown as DefaultMenuSelectOption[]}
          InputDropdownComponent={
            StyledComplexFilterInputDropdown as typeof ComplexFilterInputDropdown
          }
          DropdownMenuProps={DropdownMenuProps}
          InputDropdownProps={InputDropdownProps}
          PopperComponent={StyledPopper}
        />
        {/* (alec) disable development stage filter for now */}
        {false && (
          <StyledComplexFilter
            multiple
            data-testid="de-development-filter"
            search
            label="Development Stage"
            options={
              development_stage_terms as unknown as DefaultMenuSelectOption[]
            }
            onChange={handleDevelopmentStagesChange}
            value={
              selectedDevelopmentStages as unknown as DefaultMenuSelectOption[]
            }
            InputDropdownComponent={
              StyledComplexFilterInputDropdown as typeof ComplexFilterInputDropdown
            }
            DropdownMenuProps={DropdownMenuProps}
            InputDropdownProps={InputDropdownProps}
            PopperComponent={StyledPopper}
          />
        )}
        <StyledComplexFilter
          multiple
          data-testid="de-disease-filter"
          search
          label="Disease"
          options={disease_terms as unknown as DefaultMenuSelectOption[]}
          onChange={handleDiseasesChange}
          value={selectedDiseases as unknown as DefaultMenuSelectOption[]}
          InputDropdownComponent={
            StyledComplexFilterInputDropdown as typeof ComplexFilterInputDropdown
          }
          DropdownMenuProps={DropdownMenuProps}
          InputDropdownProps={InputDropdownProps}
          PopperComponent={StyledPopper}
        />
        <StyledComplexFilter
          multiple
          data-testid="de-self-reported-ethnicity-filter"
          search
          label="Self-Reported Ethnicity"
          options={
            self_reported_ethnicity_terms as unknown as DefaultMenuSelectOption[]
          }
          onChange={handleEthnicitiesChange}
          value={selectedEthnicities as unknown as DefaultMenuSelectOption[]}
          InputDropdownComponent={
            StyledComplexFilterInputDropdown as typeof ComplexFilterInputDropdown
          }
          DropdownMenuProps={DropdownMenuProps}
          InputDropdownProps={InputDropdownProps}
          PopperComponent={StyledPopper}
        />
        <StyledComplexFilter
          multiple
          data-testid="de-sex-filter"
          search
          label="Sex"
          options={sex_terms as unknown as DefaultMenuSelectOption[]}
          onChange={handleSexesChange}
          value={selectedSexes as unknown as DefaultMenuSelectOption[]}
          InputDropdownComponent={
            StyledComplexFilterInputDropdown as typeof ComplexFilterInputDropdown
          }
          DropdownMenuProps={DropdownMenuProps}
          InputDropdownProps={InputDropdownProps}
          PopperComponent={StyledPopper}
        />
      </div>
    </Wrapper>
  );
});

function getOptionSelected(
  option: { id: string },
  value: { id: string }
): boolean {
  return option.id === value.id;
}

function sortOptions(a: DefaultMenuSelectOption, b: DefaultMenuSelectOption) {
  if (a.name < b.name) {
    return -1;
  }
  if (a.name > b.name) {
    return 1;
  }
  return 0;
}
