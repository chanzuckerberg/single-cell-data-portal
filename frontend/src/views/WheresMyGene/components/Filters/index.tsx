import { createFilterOptions } from "@mui/material";
import {
  ComplexFilterInputDropdown,
  DefaultMenuSelectOption,
  InputDropdownProps,
  Tooltip,
} from "czifui";
import isEqual from "lodash/isEqual";
import {
  Dispatch,
  memo,
  SetStateAction,
  useCallback,
  useContext,
  useEffect,
  useMemo,
} from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import {
  FilterDimensions,
  RawDataset,
  useFilterDimensions,
} from "src/common/queries/wheresMyGene";
import { DispatchContext, StateContext } from "../../common/store";
import { selectFilters } from "../../common/store/actions";
import { Filters as IFilters } from "../../common/types";
import Organism from "./components/Organism";
import Compare from "./components/Compare";
import Sort from "./components/Sort";
import {
  StyledComplexFilter,
  StyledComplexFilterInputDropdown,
  ViewOptionsLabel,
  Wrapper,
} from "./style";
import ColorScale from "./components/ColorScale";
import { ViewOptionsWrapper } from "./components/Sort/style";

const ANALYTICS_MAPPING: {
  [key in keyof IFilters]: { eventName: EVENTS; label: string };
} = {
  datasets: {
    eventName: EVENTS.FILTER_SELECT_DATASET,
    label: "dataset_name",
  },
  diseases: {
    eventName: EVENTS.FILTER_SELECT_DISEASE,
    label: "disease",
  },
  ethnicities: {
    eventName: EVENTS.FILTER_SELECT_SELF_REPORTED_ETHNICITY,
    label: "ethnicity",
  },
  sexes: {
    eventName: EVENTS.FILTER_SELECT_SEX,
    label: "gender",
  },
};

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

export interface Props {
  isLoading: boolean;
  availableFilters: Partial<FilterDimensions>;
  setAvailableFilters: Dispatch<SetStateAction<Partial<FilterDimensions>>>;
  setIsScaled: Dispatch<SetStateAction<boolean>>;
}

export default memo(function Filters({
  isLoading,
  availableFilters,
  setAvailableFilters,
  setIsScaled,
}: Props): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const state = useContext(StateContext);

  const { selectedFilters, selectedTissues, selectedGenes } = state;

  const {
    datasets: datasetIds,
    diseases,
    ethnicities,
    sexes,
  } = selectedFilters;

  const {
    data: {
      datasets: rawDatasets,
      development_stage_terms: rawDevelopmentStages,
      disease_terms: rawDiseases,
      self_reported_ethnicity_terms: rawEthnicities,
      sex_terms: rawSexes,
    },
    isLoading: rawIsLoading,
  } = useFilterDimensions();

  const isHeatmapShown = !!selectedTissues.length && !!selectedGenes.length;

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
    };

    if (isEqual(availableFilters, newAvailableFilters)) return;

    setAvailableFilters(newAvailableFilters);
  }, [
    rawDatasets,
    rawDevelopmentStages,
    rawDiseases,
    rawEthnicities,
    rawSexes,
    rawIsLoading,
    availableFilters,
    setAvailableFilters,
  ]);

  const {
    datasets = EMPTY_ARRAY,
    disease_terms = EMPTY_ARRAY,
    self_reported_ethnicity_terms = EMPTY_ARRAY,
    sex_terms = EMPTY_ARRAY,
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

        const newlySelected = options.filter(
          (selected) => !currentOptions?.includes(selected)
        );

        // If there are newly selected filters, send an analytic event for each of them
        if (newlySelected.length) {
          newlySelected.forEach((selected) => {
            const { eventName, label } = ANALYTICS_MAPPING[key]!;
            track(eventName, {
              [label]: selected.name,
            });
          });
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

  return (
    <Wrapper>
      <div>
        <StyledComplexFilter
          multiple
          data-test-id="dataset-filter"
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
        />
        <StyledComplexFilter
          multiple
          data-test-id="disease-filter"
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
        />
        <StyledComplexFilter
          multiple
          data-test-id="self-reported-ethnicity-filter"
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
        />
        <StyledComplexFilter
          multiple
          data-test-id="sex-filter"
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
        />
      </div>

      <Organism isLoading={isLoading} />

      <Tooltip
        sdsStyle="dark"
        arrow
        placement="right"
        title={"Please select at least one tissue and gene to use this option."}
        disableHoverListener={isHeatmapShown}
        disableFocusListener={isHeatmapShown}
      >
        <div>
          <Compare areFiltersDisabled={!isHeatmapShown} />
        </div>
      </Tooltip>

      <div>
        <ViewOptionsLabel>View Options</ViewOptionsLabel>
        <ViewOptionsWrapper>
          <Tooltip
            sdsStyle="dark"
            arrow
            placement="right"
            title={
              "Please select at least one tissue and gene to use this option."
            }
            disableHoverListener={isHeatmapShown}
            disableFocusListener={isHeatmapShown}
          >
            <div>
              <Sort areFiltersDisabled={!isHeatmapShown} />
            </div>
          </Tooltip>
          <ColorScale setIsScaled={setIsScaled} />
        </ViewOptionsWrapper>
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
