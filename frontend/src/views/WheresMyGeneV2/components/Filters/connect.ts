import {
  ComplexFilterProps,
  DefaultMenuSelectOption,
  InputDropdownProps,
} from "@czi-sds/components";
import { createFilterOptions } from "@mui/base";
import {
  RawDataset,
  useFilterDimensions,
} from "src/common/queries/wheresMyGene";
import { Props } from "./types";
import { useCallback, useContext, useEffect, useMemo } from "react";
import {
  DispatchContext,
  StateContext,
} from "../../../WheresMyGene/common/store";
import { isEqual } from "lodash";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { track } from "src/common/analytics";
import { Filters as IFilters } from "src/views/WheresMyGene/common/types";
import { ANALYTICS_MAPPING } from "./constants";
import { selectFilters } from "src/views/WheresMyGene/common/store/actions";
import {
  getOptionSelected,
  isOptionEqualToValue,
  mapTermToFilterOption,
  sortOptions,
} from "./utils";

export const useConnect = ({
  availableFilters,
  setAvailableFilters,
}: {
  availableFilters: Props["availableFilters"];
  setAvailableFilters: Props["setAvailableFilters"];
}) => {
  const dispatch = useContext(DispatchContext);
  const state = useContext(StateContext);

  const { selectedFilters, selectedGenes } = state;

  const {
    datasets: datasetIds,
    diseases,
    ethnicities,
    publications,
    sexes,
    tissues,
  } = selectedFilters;

  const {
    data: {
      datasets: rawDatasets,
      development_stage_terms: rawDevelopmentStages,
      disease_terms: rawDiseases,
      self_reported_ethnicity_terms: rawEthnicities,
      publication_citations: rawPublications,
      sex_terms: rawSexes,
      tissue_terms: rawTissues,
    },
    isLoading: rawIsLoading,
  } = useFilterDimensions(2);

  const isHeatmapShown = !!selectedGenes.length;

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

    const newPublications = rawPublications.map(mapTermToFilterOption);
    newPublications.sort((a, b) => a.name.localeCompare(b.name));

    const newDevelopmentStages = rawDevelopmentStages.map(
      mapTermToFilterOption
    );
    newDevelopmentStages.sort((a, b) => a.name.localeCompare(b.name));

    const newTissues = rawTissues
      .filter((tissue) => {
        // (thuang): Product requirement to exclude "cell culture" from the list
        // https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell-data-portal/2335
        return !tissue.name.includes("(cell culture)");
      })
      .map(mapTermToFilterOption);

    newTissues.sort((a, b) => a.name.localeCompare(b.name));

    const newAvailableFilters = {
      datasets: newDatasets,
      development_stage_terms: newDevelopmentStages,
      disease_terms: newDiseases,
      self_reported_ethnicity_terms: newEthnicities,
      publication_citations: newPublications,
      sex_terms: newSexes,
      tissue_terms: newTissues,
    };

    if (isEqual(availableFilters, newAvailableFilters)) return;
    setAvailableFilters(newAvailableFilters);
  }, [
    rawDatasets,
    rawDevelopmentStages,
    rawDiseases,
    rawEthnicities,
    rawPublications,
    rawSexes,
    rawIsLoading,
    availableFilters,
    setAvailableFilters,
    rawTissues,
  ]);

  const {
    datasets = EMPTY_ARRAY,
    disease_terms = EMPTY_ARRAY,
    self_reported_ethnicity_terms = EMPTY_ARRAY,
    publication_citations = EMPTY_ARRAY,
    sex_terms = EMPTY_ARRAY,
    tissue_terms = EMPTY_ARRAY,
  } = availableFilters;

  const selectedDatasets = useMemo(() => {
    return datasets.filter((dataset) => datasetIds?.includes(dataset.id));
  }, [datasets, datasetIds]);

  const selectedDiseases = useMemo(() => {
    return disease_terms.filter((disease) => diseases?.includes(disease.id));
  }, [disease_terms, diseases]);

  const selectedEthnicities = useMemo(() => {
    return self_reported_ethnicity_terms.filter(
      (ethnicity) => ethnicities?.includes(ethnicity.id)
    );
  }, [self_reported_ethnicity_terms, ethnicities]);

  const selectedPublications = useMemo(() => {
    return publication_citations.filter(
      (publication) => publications?.includes(publication.id)
    );
  }, [publication_citations, publications]);

  const selectedSexes = useMemo(() => {
    return sex_terms.filter((sex) => sexes?.includes(sex.id));
  }, [sex_terms, sexes]);

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

        const newlySelected = options.filter(
          (selected) => !currentOptions?.includes(selected)
        );

        if (newlySelected.length) {
          newlySelected.forEach((selected) => {
            const analyticsMapping = ANALYTICS_MAPPING[key];
            if (analyticsMapping) {
              const { eventName, label } = analyticsMapping;
              track(eventName, {
                [label]: selected.name,
              });
            }
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

  const handlePublicationsChange = useMemo(
    () => handleFilterChange("publications"),
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

  const filterOptions = createFilterOptions({
    stringify: (option: RawDataset) =>
      `${option.label} ${option.collection_label}`,
  });

  const DropdownMenuProps = {
    filterOptions,
    getOptionSelected,
    isOptionEqualToValue,
  } as ComplexFilterProps<true>["DropdownMenuProps"];

  const InputDropdownProps = {
    sdsStyle: "minimal",
  } as Partial<InputDropdownProps>;

  return {
    handle: {
      tissuesChange: handleTissuesChange,
      sexesChange: handleSexesChange,
      publicationsChange: handlePublicationsChange,
      ethnicitiesChange: handleEthnicitiesChange,
      diseasesChange: handleDiseasesChange,
      datasetsChange: handleDatasetsChange,
    },
    selected: {
      tissues: selectedTissues,
      sexes: selectedSexes,
      publications: selectedPublications,
      ethnicities: selectedEthnicities,
      diseases: selectedDiseases,
      datasets: selectedDatasets,
    },
    terms: {
      tissue: tissue_terms,
      sex: sex_terms,
      publication: publication_citations,
      self_reported_ethnicity: self_reported_ethnicity_terms,
      disease: disease_terms,
    },
    DropdownMenuProps,
    isHeatmapShown,
    InputDropdownProps,
    datasets,
  };
};
