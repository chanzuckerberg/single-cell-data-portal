import { useState } from "react";
import { createFilterOptions } from "@mui/material";
import {
  ComplexFilterInputDropdown,
  DefaultMenuSelectOption,
  InputDropdownProps,
} from "@czi-sds/components";
import isEqual from "lodash/isEqual";
import { memo, useCallback, useContext, useEffect, useMemo } from "react";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import {
  FilterDimensions,
  RawDataset,
  useQueryGroupFilterDimensions,
} from "src/common/queries/differentialExpression";
import { DispatchContext } from "src/views/DifferentialExpression/common/store";
import {
  StyledComplexFilter,
  StyledComplexFilterInputDropdown,
  Wrapper,
  StyledPopper,
  StyledTagFilter,
  TagWrapper,
  EmptyRectangle,
  ClearButtonWrapper,
} from "./style";
import { QueryGroup } from "src/views/DifferentialExpression/common/store/reducer";
import {
  selectQueryGroup1Filters,
  selectQueryGroup2Filters,
  clearQueryGroup1Filters,
  clearQueryGroup2Filters,
} from "src/views/DifferentialExpression/common/store/actions";
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

interface Props {
  queryGroup: QueryGroup;
  queryGroupWithNames: QueryGroup;
  isQueryGroup1: boolean;
}
export default memo(function Filters({
  queryGroup,
  queryGroupWithNames,
  isQueryGroup1,
}: Props): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const selectQueryGroupFilters = isQueryGroup1
    ? selectQueryGroup1Filters
    : selectQueryGroup2Filters;
  const clearQueryGroupFilters = isQueryGroup1
    ? clearQueryGroup1Filters
    : clearQueryGroup2Filters;
  const [availableFilters, setAvailableFilters] =
    useState<Partial<FilterDimensions>>(EMPTY_OBJECT);

  const {
    datasets: datasetIds,
    diseases,
    ethnicities,
    sexes,
    tissues,
    cellTypes,
  } = queryGroup;

  const {
    data: {
      datasets: rawDatasets,
      development_stage_terms: rawDevelopmentStages,
      disease_terms: rawDiseases,
      self_reported_ethnicity_terms: rawEthnicities,
      sex_terms: rawSexes,
      tissue_terms: rawTissues,
      cell_type_terms: rawCellTypes,
    },
    isLoading: rawIsLoading,
  } = useQueryGroupFilterDimensions(queryGroup);

  const InputDropdownProps = {
    sdsStyle: "minimal",
  } as Partial<InputDropdownProps>;

  // (thuang): We only update available filters when API call is done,
  // otherwise when `useFilterDimensions(queryGroupIndex)` is still loading, its filters
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

    const newCellTypes = rawCellTypes.map(mapTermToFilterOption);
    newCellTypes.sort((a, b) => a.name.localeCompare(b.name));

    const newEthnicities = rawEthnicities.map(mapTermToFilterOption);
    newEthnicities.sort((a, b) => a.name.localeCompare(b.name));

    const newDevelopmentStages = rawDevelopmentStages;

    const newAvailableFilters = {
      datasets: newDatasets,
      development_stage_terms: newDevelopmentStages,
      disease_terms: newDiseases,
      self_reported_ethnicity_terms: newEthnicities,
      sex_terms: newSexes,
      tissue_terms: newTissues,
      cell_type_terms: newCellTypes,
    };

    if (isEqual(availableFilters, newAvailableFilters)) return;

    setAvailableFilters(newAvailableFilters);
  }, [
    rawDatasets,
    rawDiseases,
    rawEthnicities,
    rawSexes,
    rawCellTypes,
    rawIsLoading,
    rawDevelopmentStages,
    rawTissues,
    availableFilters,
    setAvailableFilters,
  ]);

  const {
    datasets = EMPTY_ARRAY,
    disease_terms = EMPTY_ARRAY,
    self_reported_ethnicity_terms = EMPTY_ARRAY,
    sex_terms = EMPTY_ARRAY,
    tissue_terms = EMPTY_ARRAY,
    cell_type_terms = EMPTY_ARRAY,
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

  const selectedTissues = useMemo(() => {
    return tissue_terms.filter((tissue) => tissues?.includes(tissue.id));
  }, [tissue_terms, tissues]);

  const selectedCellTypes = useMemo(() => {
    return cell_type_terms.filter((cellType) =>
      cellTypes?.includes(cellType.id)
    );
  }, [cell_type_terms, cellTypes]);

  const handleFilterChange = useCallback(
    function handleFilterChange_(
      key: keyof QueryGroup
    ): (options: DefaultMenuSelectOption[] | null) => void {
      let currentOptions: DefaultMenuSelectOption[] | null = null;

      return (options: DefaultMenuSelectOption[] | null): void => {
        if (
          !dispatch ||
          !options ||
          !rawIsLoading ||
          // If the options are the same
          JSON.stringify(options.sort(sortOptions)) ===
            JSON.stringify(currentOptions?.sort(sortOptions)) ||
          // If the options change from null to [], which is the default value
          (currentOptions === null && JSON.stringify(options) === "[]")
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
    [dispatch, selectQueryGroupFilters, rawIsLoading]
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

  const handleCellTypesChange = useMemo(
    () => handleFilterChange("cellTypes"),
    [handleFilterChange]
  );

  const tagsToShow = [];
  const deleteHandlers: (() => void)[] = [];
  for (const key in queryGroupWithNames) {
    for (const [index, value] of queryGroupWithNames[
      key as keyof QueryGroup
    ].entries()) {
      tagsToShow.push(value);
      deleteHandlers.push(() => {
        if (!dispatch) return;
        const newOptionsWithNames = queryGroupWithNames[
          key as keyof QueryGroup
        ].filter((_, i) => i !== index);
        const newOptions = queryGroup[key as keyof QueryGroup].filter(
          (_, i) => i !== index
        );
        const options = newOptions.map((option, i) => {
          return { id: option, name: newOptionsWithNames[i] };
        });
        dispatch(selectQueryGroupFilters(key as keyof QueryGroup, options));
      });
    }
  }

  const handleClearQueryGroup = () => {
    if (!dispatch) return;
    dispatch(clearQueryGroupFilters());
  };

  const isActive = !!tagsToShow.length;
  return (
    <Wrapper>
      <ClearButtonWrapper onClick={handleClearQueryGroup}>
        Clear
      </ClearButtonWrapper>
      <TagWrapper>
        {isActive ? (
          tagsToShow.map((tag, index) => (
            <StyledTagFilter
              key={index}
              onDelete={deleteHandlers[index]}
              label={tag}
            />
          ))
        ) : (
          <EmptyRectangle />
        )}
      </TagWrapper>
      <StyledComplexFilter
        multiple
        data-testid="de-qg-tissue-filter"
        search
        label={"Tissue"}
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
        data-testid="de-qg-cell-type-filter"
        search
        label="Cell Type"
        options={cell_type_terms as unknown as DefaultMenuSelectOption[]}
        onChange={handleCellTypesChange}
        value={selectedCellTypes as unknown as DefaultMenuSelectOption[]}
        InputDropdownComponent={
          StyledComplexFilterInputDropdown as typeof ComplexFilterInputDropdown
        }
        DropdownMenuProps={DropdownMenuProps}
        InputDropdownProps={InputDropdownProps}
        PopperComponent={StyledPopper}
      />
      <StyledComplexFilter
        multiple
        data-testid="de-qg-dataset-filter"
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
      <StyledComplexFilter
        multiple
        data-testid="de-qg-disease-filter"
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
        data-testid="de-qg-self-reported-ethnicity-filter"
        search
        label="Ethnicity"
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
        data-testid="de-qg-sex-filter"
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
