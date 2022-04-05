import {
  ComplexFilter,
  ComplexFilterInputDropdown,
  DefaultMenuSelectOption,
} from "czifui";
import isEqual from "lodash/isEqual";
import {
  memo,
  useCallback,
  useContext,
  useEffect,
  useMemo,
  useState,
} from "react";
import { EMPTY_ARRAY, EMPTY_OBJECT } from "src/common/constants/utils";
import {
  FilterDimensions,
  useFilterDimensions,
} from "src/common/queries/wheresMyGene";
import { DispatchContext, StateContext } from "../../common/store";
import { selectFilters } from "../../common/store/actions";
import { Filters as IFilters } from "../../common/types";
import { StyledComplexFilterInputDropdown } from "./style";

const MenuSelectProps = {
  getOptionSelected,
};

export default memo(function Filters(): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const state = useContext(StateContext);
  const [availableFilters, setAvailableFilters] =
    useState<Partial<FilterDimensions>>(EMPTY_OBJECT);

  const { selectedFilters } = state;

  const {
    datasets: datasetIds,
    developmentStages,
    diseases,
    ethnicities,
    sexes,
  } = selectedFilters;

  const {
    data: {
      datasets: rawDatasets,
      development_stage_terms: rawDevelopmentStages,
      disease_terms: rawDiseases,
      ethnicity_terms: rawEthnicities,
      sex_terms: rawSexes,
    },
    isLoading: rawIsLoading,
  } = useFilterDimensions({ includeAllFilterOptions: true });

  // (thuang): We only update available filters when API call is done,
  // otherwise when `useFilterDimensions()` is still loading, its filters
  // will temporarily be empty, and thus resetting the selected filter values
  useEffect(() => {
    if (rawIsLoading) return;

    const newAvailableFilters = {
      datasets: rawDatasets.map((dataset) => ({
        ...dataset,
        name: dataset.label,
      })),
      development_stage_terms: rawDevelopmentStages,
      disease_terms: rawDiseases,
      ethnicity_terms: rawEthnicities,
      sex_terms: rawSexes,
    };

    if (isEqual(availableFilters, newAvailableFilters)) return;

    setAvailableFilters(newAvailableFilters);
  }, [rawDatasets, rawDevelopmentStages, rawDiseases, rawEthnicities, rawSexes, rawIsLoading, availableFilters, setAvailableFilters]);

  const {
    datasets = EMPTY_ARRAY,
    development_stage_terms = EMPTY_ARRAY,
    disease_terms = EMPTY_ARRAY,
    ethnicity_terms = EMPTY_ARRAY,
    sex_terms = EMPTY_ARRAY,
  } = availableFilters;

  const selectedDatasets = useMemo(() => {
    return datasets.filter((dataset) => datasetIds?.includes(dataset.id));
  }, [datasets, datasetIds]);

  const selectedDevelopmentStages = useMemo(() => {
    return development_stage_terms.filter((developmentStage) =>
      developmentStages?.includes(developmentStage.id)
    );
  }, [development_stage_terms, developmentStages]);

  const selectedDiseases = useMemo(() => {
    return disease_terms.filter((disease) => diseases?.includes(disease.id));
  }, [disease_terms, diseases]);

  const selectedEthnicities = useMemo(() => {
    return ethnicity_terms.filter((ethnicity) =>
      ethnicities?.includes(ethnicity.id)
    );
  }, [ethnicity_terms, ethnicities]);

  const selectedSexes = useMemo(() => {
    return sex_terms.filter((sex) => sexes?.includes(sex.id));
  }, [sex_terms, sexes]);

  const handleFilterChange = useCallback(
    function handleFilterChange(
      key: keyof IFilters
    ): (options: DefaultMenuSelectOption[] | null) => void {
      return (options: DefaultMenuSelectOption[] | null): void => {
        if (!dispatch || !options) return;

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

  const handleDevelopmentStagesChange = useMemo(
    () => handleFilterChange("developmentStages"),
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

  return (
    <div>
      <ComplexFilter
        multiple
        label="Dataset"
        options={datasets as unknown as DefaultMenuSelectOption[]}
        onChange={handleDatasetsChange}
        value={selectedDatasets as unknown as DefaultMenuSelectOption[]}
        InputDropdownComponent={
          StyledComplexFilterInputDropdown as typeof ComplexFilterInputDropdown
        }
        MenuSelectProps={MenuSelectProps}
      />
      <ComplexFilter
        multiple
        label="Development Stage"
        options={development_stage_terms}
        onChange={handleDevelopmentStagesChange}
        value={selectedDevelopmentStages}
        InputDropdownComponent={
          StyledComplexFilterInputDropdown as typeof ComplexFilterInputDropdown
        }
        MenuSelectProps={MenuSelectProps}
      />
      <ComplexFilter
        multiple
        label="Disease"
        options={disease_terms}
        onChange={handleDiseasesChange}
        value={selectedDiseases}
        InputDropdownComponent={
          StyledComplexFilterInputDropdown as typeof ComplexFilterInputDropdown
        }
        MenuSelectProps={MenuSelectProps}
      />
      <ComplexFilter
        multiple
        label="Ethnicity"
        options={ethnicity_terms}
        onChange={handleEthnicitiesChange}
        value={selectedEthnicities}
        InputDropdownComponent={
          StyledComplexFilterInputDropdown as typeof ComplexFilterInputDropdown
        }
        MenuSelectProps={MenuSelectProps}
      />
      <ComplexFilter
        multiple
        label="Sex"
        options={sex_terms}
        onChange={handleSexesChange}
        value={selectedSexes}
        InputDropdownComponent={
          StyledComplexFilterInputDropdown as typeof ComplexFilterInputDropdown
        }
        MenuSelectProps={MenuSelectProps}
      />
    </div>
  );
});

function getOptionSelected(
  option: { id: string },
  value: { id: string }
): boolean {
  return option.id === value.id;
}
