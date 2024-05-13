import { useEffect, useState } from "react";
import isEqual from "lodash/isEqual";
import { useQueryGroupFilterDimensions } from "src/common/queries/differentialExpression";
import { QueryGroup } from "src/views/DifferentialExpression/common/store/reducer";
import { FilterOption } from "../Filters/types";

export const EMPTY_FILTERS = {
  disease_terms: [],
  self_reported_ethnicity_terms: [],
  sex_terms: [],
  tissue_terms: [],
  cell_type_terms: [],
  publication_citations: [],
};

export interface FilterOptionDimensions {
  disease_terms: FilterOption[];
  self_reported_ethnicity_terms: FilterOption[];
  sex_terms: FilterOption[];
  tissue_terms: FilterOption[];
  cell_type_terms: FilterOption[];
  publication_citations: FilterOption[];
}

const _mapTermToFilterOption = (term: {
  id: string;
  name: string;
}): FilterOption => {
  return {
    name: term.name,
    id: term.id,
  };
};

function useProcessedQueryGroupFilterDimensions(queryGroup: QueryGroup): {
  availableFilters: FilterOptionDimensions;
  n_cells: number;
} {
  const [localNCells, setLocalNCells] = useState<number>(0);

  const [availableFilters, setAvailableFilters] =
    useState<FilterOptionDimensions>(EMPTY_FILTERS);

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
    n_cells,
    isLoading: rawIsLoading,
  } = useQueryGroupFilterDimensions(queryGroup);

  useEffect(() => {
    if (rawIsLoading) return;
    const newPublicationCitations = rawPublicationCitations.map((citation) => ({
      name: citation,
      id: citation,
    }));
    setLocalNCells(n_cells);

    newPublicationCitations.sort((a, b) => a.name.localeCompare(b.name));

    const newSexes = rawSexes.map(_mapTermToFilterOption);
    newSexes.sort((a, b) => a.name.localeCompare(b.name));

    const newDiseases = rawDiseases.map(_mapTermToFilterOption);
    newDiseases.sort((a, b) =>
      a.name === "normal"
        ? -1
        : b.name === "normal"
        ? 1
        : a.name.localeCompare(b.name)
    );

    const newTissues = rawTissues.map(_mapTermToFilterOption);
    newTissues.sort((a, b) => a.name.localeCompare(b.name));
    console.log(newTissues);
    const newCellTypes = rawCellTypes.map(_mapTermToFilterOption);
    newCellTypes.sort((a, b) => a.name.localeCompare(b.name));

    const newEthnicities = rawEthnicities.map(_mapTermToFilterOption);
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
    n_cells,
  ]);

  return { availableFilters, n_cells: localNCells };
}

export default useProcessedQueryGroupFilterDimensions;
