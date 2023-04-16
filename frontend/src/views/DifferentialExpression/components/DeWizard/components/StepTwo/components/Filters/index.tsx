import { useState } from "react";
import { createFilterOptions } from "@mui/material";
import {
  ComplexFilterInputDropdown,
  DefaultMenuSelectOption,
  Icon,
  InputDropdownProps,
} from "czifui";
import isEqual from "lodash/isEqual";
import { memo, useCallback, useContext, useEffect, useMemo } from "react";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import {
  FilterDimensions,
  RawDataset,
  useQueryGroupFilterDimensions,
} from "src/common/queries/differentialExpression";
import {
  DispatchContext,
  StateContext,
} from "src/views/DifferentialExpression/common/store";
import {
  selectQueryGroupFilters,
  deleteQueryGroup,
} from "src/views/DifferentialExpression/common/store/actions";
import {
  StyledComplexFilter,
  StyledComplexFilterInputDropdown,
  Wrapper,
  StyledPopper,
  StyledTagFilter,
  TagWrapper,
  QueryGroupTitle,
  EmptyRectangle,
  IconButtonWrapper,
} from "./style";
import {
  QueryGroup,
  QueryGroupWithNames,
} from "src/views/DifferentialExpression/common/store/reducer";

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
  queryGroupIndex: number;
  queryGroup: QueryGroup;
  queryGroupWithNames: QueryGroupWithNames;
}
export default memo(function Filters({
  queryGroupIndex,
  queryGroup,
  queryGroupWithNames,
}: Props): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const state = useContext(StateContext);
  const { selectedFilters } = state;

  const [availableFilters, setAvailableFilters] =
    useState<Partial<FilterDimensions>>(EMPTY_OBJECT);

  const {
    datasets: datasetIds,
    diseases,
    ethnicities,
    sexes,
    tissues,
    developmentStages,
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
  } = useQueryGroupFilterDimensions(queryGroup, availableFilters);

  const InputDropdownProps = {
    sdsStyle: "minimal",
  } as Partial<InputDropdownProps>;

  // (thuang): We only update available filters when API call is done,
  // otherwise when `useFilterDimensions(queryGroupIndex)` is still loading, its filters
  // will temporarily be empty, and thus resetting the selected filter values
  useEffect(() => {
    if (rawIsLoading) return;
    const newDatasets = rawDatasets
      .filter((dataset) => {
        return selectedFilters.datasets.length
          ? selectedFilters.datasets.includes(dataset.id)
          : true;
      })
      .map((dataset) => ({
        ...dataset,
        details: dataset.collection_label,
        name: dataset.label,
      }));
    newDatasets.sort((a, b) => a.name.localeCompare(b.name));

    const newSexes = rawSexes
      .filter((sex) => {
        return selectedFilters.sexes.length
          ? selectedFilters.sexes.includes(sex.id)
          : true;
      })
      .map(mapTermToFilterOption);
    newSexes.sort((a, b) => a.name.localeCompare(b.name));

    const newDiseases = rawDiseases
      .filter((disease) => {
        return selectedFilters.diseases.length
          ? selectedFilters.diseases.includes(disease.id)
          : true;
      })
      .map(mapTermToFilterOption);
    newDiseases.sort((a, b) =>
      a.name === "normal"
        ? -1
        : b.name === "normal"
        ? 1
        : a.name.localeCompare(b.name)
    );

    const newTissues = rawTissues
      .filter((tissue) => {
        return selectedFilters.tissues.length
          ? selectedFilters.tissues.includes(tissue.id)
          : true;
      })
      .map(mapTermToFilterOption);
    newTissues.sort((a, b) => a.name.localeCompare(b.name));

    const newCellTypes = rawCellTypes.map(mapTermToFilterOption);
    newCellTypes.sort((a, b) => a.name.localeCompare(b.name));

    const newEthnicities = rawEthnicities
      .filter((ethnicity) => {
        return selectedFilters.ethnicities.length
          ? selectedFilters.ethnicities.includes(ethnicity.id)
          : true;
      })
      .map(mapTermToFilterOption);
    newEthnicities.sort((a, b) => a.name.localeCompare(b.name));

    const newDevelopmentStages = rawDevelopmentStages
      .filter((stage) => {
        return selectedFilters.developmentStages.length
          ? selectedFilters.developmentStages.includes(stage.id)
          : true;
      })
      .map(mapTermToFilterOption);
    newDevelopmentStages.sort((a, b) => a.name.localeCompare(b.name));

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
    rawDevelopmentStages,
    rawDiseases,
    rawEthnicities,
    rawSexes,
    rawCellTypes,
    rawIsLoading,
    availableFilters,
    setAvailableFilters,
  ]);

  const {
    datasets = EMPTY_ARRAY,
    disease_terms = EMPTY_ARRAY,
    self_reported_ethnicity_terms = EMPTY_ARRAY,
    sex_terms = EMPTY_ARRAY,
    tissue_terms = EMPTY_ARRAY,
    development_stage_terms = EMPTY_ARRAY,
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

  const selectedDevelopmentStages = useMemo(() => {
    return development_stage_terms.filter((stage) =>
      developmentStages?.includes(stage.id)
    );
  }, [development_stage_terms, developmentStages]);

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
        dispatch(
          selectQueryGroupFilters(key, optionsWithNames, queryGroupIndex)
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

  const handleCellTypesChange = useMemo(
    () => handleFilterChange("cellTypes"),
    [handleFilterChange]
  );

  const tagsToShow = [];
  const deleteHandlers: (() => void)[] = [];
  for (const key in queryGroupWithNames) {
    for (const [index, value] of queryGroupWithNames[
      key as keyof QueryGroupWithNames
    ].entries()) {
      tagsToShow.push(value);
      deleteHandlers.push(() => {
        if (!dispatch) return;
        const newOptionsWithNames = queryGroupWithNames[
          key as keyof QueryGroupWithNames
        ].filter((_, i) => i !== index);
        const newOptions = queryGroup[key as keyof QueryGroup].filter(
          (_, i) => i !== index
        );
        const options = newOptions.map((option, i) => {
          return { id: option, name: newOptionsWithNames[i] };
        });
        dispatch(
          selectQueryGroupFilters(
            key as keyof QueryGroup,
            options,
            queryGroupIndex
          )
        );
      });
    }
  }

  const handleDeleteQueryGroup = () => {
    if (!dispatch) return;
    dispatch(deleteQueryGroup(queryGroupIndex));
  };
  const isActive = !!tagsToShow.length;
  return (
    <Wrapper active={isActive}>
      <QueryGroupTitle>
        Query group
        <IconButtonWrapper onClick={handleDeleteQueryGroup}>
          <Icon sdsIcon="trashCan" sdsSize="s" sdsType="button" />
        </IconButtonWrapper>
      </QueryGroupTitle>
      <TagWrapper>
        {isActive ? (
          tagsToShow.map((tag, index) => (
            <StyledTagFilter onDelete={deleteHandlers[index]} label={tag} />
          ))
        ) : (
          <EmptyRectangle />
        )}
      </TagWrapper>

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
        data-testid="de-qg-tissue-filter"
        search
        label={"Tissue\u002A"}
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
      {/* (alec) disable development stage filter for now */}
      {false && (
        <StyledComplexFilter
          multiple
          data-testid="de-qg-development-filter"
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
