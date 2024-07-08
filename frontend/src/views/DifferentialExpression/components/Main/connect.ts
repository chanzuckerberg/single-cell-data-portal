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
import { craftPayloadWithQueryGroups } from "./utils";
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

  const handleRunDifferentialExpression = () => {
    if (!dispatch) return;
    dispatch(submitQueryGroups());

    track(
      EVENTS.DE_FIND_GENES_CLICKED,
      craftPayloadWithQueryGroups(queryGroups)
    );
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
    if (!organism) return;

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
    const cellTypes = params.get("celltypes")?.split(",");
    const diseases = params.get("diseases")?.split(",");
    const ethnicities = params.get("ethnicities")?.split(",");
    const publications = params.get("publications")?.split(",");
    const sexes = params.get("sexes")?.split(",");
    const tissues = params.get("tissues")?.split(",");
    const organism: string | null = params.get("organism");

    if (organism) {
      paramsToRemove.push("organism");
    }

    if (tissues) {
      paramsToRemove.push("tissues");
      const tissuesFiltered = allFilterOptions.tissue_terms.filter((tissue) =>
        tissues.includes(tissue.id)
      );
      dispatch(selectQueryGroup1Filters("tissues", tissuesFiltered));
    }

    if (cellTypes) {
      paramsToRemove.push("celltypes");
      const cellTypesFiltered = allFilterOptions.cell_type_terms.filter(
        (cellType) => cellTypes.includes(cellType.id)
      );
      dispatch(selectQueryGroup1Filters("cellTypes", cellTypesFiltered));
    }
    if (diseases) {
      paramsToRemove.push("diseases");
      const diseasesFiltered = allFilterOptions.disease_terms.filter(
        (disease) => diseases.includes(disease.id)
      );
      dispatch(selectQueryGroup1Filters("diseases", diseasesFiltered));
    }
    if (ethnicities) {
      paramsToRemove.push("ethnicities");
      const ethnicitiesFiltered =
        allFilterOptions.self_reported_ethnicity_terms.filter((ethnicity) =>
          ethnicities.includes(ethnicity.id)
        );
      dispatch(selectQueryGroup1Filters("ethnicities", ethnicitiesFiltered));
    }
    if (publications) {
      paramsToRemove.push("publications");
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
      paramsToRemove.push("sexes");
      const sexesFiltered = allFilterOptions.sex_terms.filter((sex) =>
        sexes.includes(sex.id)
      );
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
