import { useContext, useEffect, useState } from "react";
import {
  DispatchContext,
  StateContext,
} from "src/views/DifferentialExpression/common/store";

import {
  submitQueryGroups,
  clearQueryGroup1Filters,
  clearQueryGroup2Filters,
  selectQueryGroup1Filters,
  selectOrganism,
} from "src/views/DifferentialExpression/common/store/actions";
import useProcessedQueryGroupFilterDimensions from "./components/common/query_group_filter_dimensions";

import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { useCraftPayloadWithQueryGroups } from "./utils";
import { removeParams } from "src/common/utils/removeParams";
import { useRouter } from "next/router";
import {
  useAvailableOrganisms,
  useQueryGroupFilterDimensions,
} from "src/common/queries/differentialExpression";
import { EMPTY_FILTERS } from "../../common/store/reducer";

export const useConnect = () => {
  // These flags control a state machine that ensures proper initialization of the DE state
  // from the URL query parameters. First organism is initialized, then the query group is initialized.
  const [organismInitializedFromShareURL, setOrganismInitializedFromShareURL] =
    useState<boolean>(false);
  const [
    queryGroupInitializedFromShareURL,
    setQueryGroupInitializedFromShareURL,
  ] = useState<boolean>(false);

  const [isLoading, setIsLoading] = useState<boolean>(false);
  const [isLoadingGetDeQuery, setIsLoadingGetDeQuery] =
    useState<boolean>(false);

  useEffect(() => {
    setIsLoadingGetDeQuery(isLoadingGetDeQuery);
  }, [isLoadingGetDeQuery]);
  const dispatch = useContext(DispatchContext);
  const { queryGroups } = useContext(StateContext);
  const { queryGroup1, queryGroup2 } = queryGroups;

  // check if any values in queryGroup1 are not empty
  const isQueryGroup1NotEmpty = Object.values(queryGroup1).some(
    (value) => value.length > 0
  );
  const isQueryGroup2NotEmpty = Object.values(queryGroup2).some(
    (value) => value.length > 0
  );
  const canRunDifferentialExpression =
    !isLoading && isQueryGroup1NotEmpty && isQueryGroup2NotEmpty;

  const payload = useCraftPayloadWithQueryGroups();

  const handleRunDifferentialExpression = () => {
    if (!dispatch) return;
    dispatch(submitQueryGroups());

    track(EVENTS.DE_FIND_GENES_CLICKED, payload);
  };

  const handleClearQueryGroups = () => {
    if (!dispatch) return;
    dispatch(clearQueryGroup1Filters());
    dispatch(clearQueryGroup2Filters());
  };

  const { n_cells: nCellsGroup1, isLoading: isLoadingGroup1 } =
    useProcessedQueryGroupFilterDimensions(queryGroup1);
  const { n_cells: nCellsGroup2, isLoading: isLoadingGroup2 } =
    useProcessedQueryGroupFilterDimensions(queryGroup2);

  const router = useRouter();
  const { data: availableOrganisms, isLoading: isLoadingOrganisms } =
    useAvailableOrganisms();

  const { data: allFilterOptions, isLoading: isLoadingFilterOptions } =
    useQueryGroupFilterDimensions(EMPTY_FILTERS);

  // We have a state machine that ensures the below effects are only triggered on page load.
  // First select the organism and set organismInitializedFromShareURL to true.
  useEffect(() => {
    if (!dispatch || isLoadingOrganisms || organismInitializedFromShareURL)
      return;
    const { search } = window.location;
    const params = new URLSearchParams(search);
    const organism: string | null = params.get("organism");

    if (!organism) {
      setOrganismInitializedFromShareURL(true); // homo sapiens is default
      return;
    }

    const isOrganismValid = availableOrganisms.some(
      (org) => org.id === organism
    );
    if (isOrganismValid) {
      dispatch(selectOrganism(organism));
    }
    setOrganismInitializedFromShareURL(true);
  }, [
    dispatch,
    isLoadingOrganisms,
    availableOrganisms,
    organismInitializedFromShareURL,
  ]);

  // Now, we parse the query groups from the URL and set queryGroupInitializedFromShareURL to true.
  // We only want to do this if the organism has been initialized from the URL since the valid filter options
  // will change based on the organism.
  useEffect(() => {
    if (
      !dispatch ||
      !organismInitializedFromShareURL ||
      isLoadingFilterOptions ||
      queryGroupInitializedFromShareURL
    )
      return;

    const { search } = window.location;
    const params = new URLSearchParams(search);

    const paramsToRemove: string[] = [];

    const getParamsAndRemove = (paramName: string) => {
      const paramValues = params.get(paramName)?.split(",");
      if (paramValues) {
        paramsToRemove.push(paramName);
      }
      return paramValues;
    };

    const filterOptions = (
      paramValues: string[],
      filterTerms: { id: string; name: string }[]
    ) => {
      return filterTerms.filter((term) => paramValues.includes(term.id));
    };

    const cellTypes = getParamsAndRemove("celltypes");
    const diseases = getParamsAndRemove("diseases");
    const ethnicities = getParamsAndRemove("ethnicities");
    const publications = getParamsAndRemove("publications");
    const sexes = getParamsAndRemove("sexes");
    const tissues = getParamsAndRemove("tissues");
    const organism: string | null = params.get("organism");

    if (organism) {
      paramsToRemove.push("organism");
    }

    if (tissues) {
      const tissuesFiltered = filterOptions(
        tissues,
        allFilterOptions.tissue_terms
      );
      dispatch(selectQueryGroup1Filters("tissues", tissuesFiltered));
    }

    if (cellTypes) {
      const cellTypesFiltered = filterOptions(
        cellTypes,
        allFilterOptions.cell_type_terms
      );
      dispatch(selectQueryGroup1Filters("cellTypes", cellTypesFiltered));
    }
    if (diseases) {
      const diseasesFiltered = filterOptions(
        diseases,
        allFilterOptions.disease_terms
      );
      dispatch(selectQueryGroup1Filters("diseases", diseasesFiltered));
    }
    if (ethnicities) {
      const ethnicitiesFiltered = filterOptions(
        ethnicities,
        allFilterOptions.self_reported_ethnicity_terms
      );
      dispatch(selectQueryGroup1Filters("ethnicities", ethnicitiesFiltered));
    }

    if (publications) {
      const publicationsFiltered = allFilterOptions.publication_citations
        .filter((publication) => publications.includes(publication))
        .map((publication) => ({
          id: publication,
          name: publication,
        }));
      dispatch(
        selectQueryGroup1Filters("publicationCitations", publicationsFiltered)
      );
    }
    if (sexes) {
      const sexesFiltered = filterOptions(sexes, allFilterOptions.sex_terms);
      dispatch(selectQueryGroup1Filters("sexes", sexesFiltered));
    }
    removeParams({
      params: paramsToRemove,
      router: router,
    });
    setQueryGroupInitializedFromShareURL(true);
  }, [
    router,
    allFilterOptions,
    isLoadingFilterOptions,
    dispatch,
    organismInitializedFromShareURL,
    setOrganismInitializedFromShareURL,
    queryGroupInitializedFromShareURL,
    setQueryGroupInitializedFromShareURL,
  ]);

  return {
    isLoading,
    setIsLoading,
    queryGroup1,
    queryGroup2,
    canRunDifferentialExpression,
    handleRunDifferentialExpression,
    handleClearQueryGroups,
    nCellsGroup1,
    isLoadingGroup1,
    nCellsGroup2,
    isLoadingGroup2,
  };
};
