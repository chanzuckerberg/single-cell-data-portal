import { useContext, useMemo } from "react";
import { useQuery, UseQueryResult } from "react-query";
import { API_URL } from "src/configs/configs";
import {
  DispatchContext,
  StateContext,
} from "src/views/DifferentialExpression/common/store";
import { setSnapshotId } from "src/views/DifferentialExpression/common/store/actions";
import { API } from "../API";
import { DEFAULT_FETCH_OPTIONS, JSON_BODY_FETCH_OPTIONS } from "./common";
import { ENTITIES } from "./entities";
import {
  EMPTY_FILTERS,
  QueryGroup,
} from "src/views/DifferentialExpression/common/store/reducer";
import { RawDataset } from "./wheresMyGene";

interface RawOntologyTerm {
  [id: string]: string;
}

interface FilterSecondary {
  organism_ontology_term_id: string;
  tissue_ontology_term_ids: string[];
  disease_ontology_term_ids: string[];
  sex_ontology_term_ids: string[];
  development_stage_ontology_term_ids: string[];
  self_reported_ethnicity_ontology_term_ids: string[];
  cell_type_ontology_term_ids: string[];
  publication_citations: string[];
}

interface FiltersQuery {
  filter: FilterSecondary;
}

interface DifferentialExpressionQuery {
  queryGroup1Filters: FilterSecondary;
  queryGroup2Filters: FilterSecondary;
}

interface FiltersQueryResponse {
  filter_dims: {
    datasets: RawDataset[];
    disease_terms: { [id: string]: string }[];
    sex_terms: { [id: string]: string }[];
    development_stage_terms: { [id: string]: string }[];
    self_reported_ethnicity_terms: { [id: string]: string }[];
    tissue_terms: { [id: string]: string }[];
    cell_type_terms: { [id: string]: string }[];
    organism_terms: { [id: string]: string }[];
    publication_citations: string[];
  };
  n_cells: number;
  snapshot_id: string;
}

interface DifferentialExpressionResult {
  gene_ontology_term_id: string;
  gene_symbol: string;
  log_fold_change: number;
  effect_size: number;
  adjusted_p_value: number;
}

interface DifferentialExpressionQueryResponse {
  differentialExpressionResults: DifferentialExpressionResult[];
  successCode: number;
  snapshot_id: string;
  n_overlap: number;
}

interface DifferentialExpressionQueryResult {
  differentialExpressionResults: DifferentialExpressionResult[];
  nOverlap: DifferentialExpressionQueryResponse["n_overlap"];
}

async function fetchFiltersQuery({
  query,
  signal,
}: {
  query: FiltersQuery | null;
  signal?: AbortSignal;
}): Promise<FiltersQueryResponse | undefined> {
  if (!query) return;

  const url = API_URL + API.DE_FILTERS_QUERY;

  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    body: JSON.stringify(query),
    method: "POST",
    signal,
  });
  const json: FiltersQueryResponse = await response.json();

  if (!response.ok) {
    throw json;
  }

  return json;
}

const USE_FILTERS_QUERY = {
  entities: [ENTITIES.DE_FILTERS_QUERY],
  id: "de-filters-query",
};

async function fetchDifferentialExpressionQuery({
  query,
  signal,
}: {
  query: DifferentialExpressionQuery | null;
  signal?: AbortSignal;
}): Promise<DifferentialExpressionQueryResponse | undefined> {
  if (!query) return;

  const url = API_URL + API.DE_QUERY;

  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    body: JSON.stringify(query),
    method: "POST",
    signal,
  });
  const json: DifferentialExpressionQueryResponse = await response.json();

  if (!response.ok) {
    throw json;
  }

  return json;
}

const USE_DE_QUERY = {
  entities: [ENTITIES.DE_QUERY],
  id: "de-query",
};

function useDEFiltersQuery(
  query: FiltersQuery | null
): UseQueryResult<FiltersQueryResponse> {
  const dispatch = useContext(DispatchContext);

  // (thuang): Refresh query when the snapshotId changes
  const currentSnapshotId = useSnapshotId();

  return useQuery(
    [USE_FILTERS_QUERY, query, currentSnapshotId],
    ({ signal }) => fetchFiltersQuery({ query, signal }),
    {
      enabled: Boolean(query),
      onSuccess(response) {
        if (!response || !dispatch) return;

        const { snapshot_id } = response;

        if (currentSnapshotId !== snapshot_id) {
          dispatch(setSnapshotId(snapshot_id));
        }
      },
      // (thuang): We don't need to refetch during the session
      staleTime: Infinity,
    }
  );
}

function useDEQuery(
  query: DifferentialExpressionQuery | null
): UseQueryResult<DifferentialExpressionQueryResponse> {
  const dispatch = useContext(DispatchContext);

  // (thuang): Refresh query when the snapshotId changes
  const currentSnapshotId = useSnapshotId();

  return useQuery(
    [USE_DE_QUERY, query, currentSnapshotId],
    ({ signal }) => fetchDifferentialExpressionQuery({ query, signal }),
    {
      enabled: Boolean(query),
      onSuccess(response) {
        if (!response || !dispatch) return;

        const { snapshot_id } = response;

        if (currentSnapshotId !== snapshot_id) {
          dispatch(setSnapshotId(snapshot_id));
        }
      },
      // (thuang): We don't need to refetch during the session
      staleTime: Infinity,
    }
  );
}

const EMPTY_FILTER_DIMENSIONS = {
  datasets: [],
  publication_citations: [],
  development_stage_terms: [],
  disease_terms: [],
  self_reported_ethnicity_terms: [],
  sex_terms: [],
  tissue_terms: [],
  cell_type_terms: [],
  organism_terms: [],
};

interface FilterDimensions {
  datasets: RawDataset[];
  development_stage_terms: { id: string; name: string }[];
  disease_terms: { id: string; name: string }[];
  self_reported_ethnicity_terms: { id: string; name: string }[];
  sex_terms: { id: string; name: string }[];
  tissue_terms: { id: string; name: string }[];
  cell_type_terms: { id: string; name: string }[];
  organism_terms: { id: string; name: string }[];
  publication_citations: string[];
}

const CXG_CENSUS_ORGANISMS = [
  { name: "Homo sapiens", id: "NCBITaxon:9606" },
  { name: "Mus musculus", id: "NCBITaxon:10090" },
];

export function useAvailableOrganisms() {
  const { data, isLoading } = useQueryGroupFilterDimensions(EMPTY_FILTERS);
  return useMemo(() => {
    if (isLoading || !data) return { data: [], isLoading };

    const organism_terms = CXG_CENSUS_ORGANISMS;

    return {
      isLoading,
      data: organism_terms.filter((organism) =>
        CXG_CENSUS_ORGANISMS.map((obj) => obj.name).includes(organism.name)
      ),
    };
  }, [data, isLoading]);
}

export function useDifferentialExpression(): {
  data: DifferentialExpressionQueryResult;
  isLoading: boolean;
} {
  const requestBody = useDEQueryRequestBody();
  const { data, isLoading } = useDEQuery(requestBody);
  return useMemo(() => {
    if (isLoading || !data)
      return {
        data: {
          differentialExpressionResults: [],
          nOverlap: 0,
        },
        isLoading,
      };
    return {
      data: {
        differentialExpressionResults: data.differentialExpressionResults,
        nOverlap: data.n_overlap,
      },
      isLoading: false,
    };
  }, [data, isLoading]);
}

function useDEQueryRequestBody() {
  const { organismId, submittedQueryGroups: queryGroups } =
    useContext(StateContext);

  return useMemo(() => {
    if (!queryGroups || !organismId) {
      return null;
    }

    const { queryGroup1, queryGroup2 } = queryGroups;

    return {
      queryGroup1Filters: {
        development_stage_ontology_term_ids: queryGroup1.developmentStages,
        disease_ontology_term_ids: queryGroup1.diseases,
        organism_ontology_term_id: organismId,
        self_reported_ethnicity_ontology_term_ids: queryGroup1.ethnicities,
        sex_ontology_term_ids: queryGroup1.sexes,
        tissue_ontology_term_ids: queryGroup1.tissues,
        cell_type_ontology_term_ids: queryGroup1.cellTypes,
        publication_citations: queryGroup1.publicationCitations,
      },
      queryGroup2Filters: {
        development_stage_ontology_term_ids: queryGroup2.developmentStages,
        disease_ontology_term_ids: queryGroup2.diseases,
        organism_ontology_term_id: organismId,
        self_reported_ethnicity_ontology_term_ids: queryGroup2.ethnicities,
        sex_ontology_term_ids: queryGroup2.sexes,
        tissue_ontology_term_ids: queryGroup2.tissues,
        cell_type_ontology_term_ids: queryGroup2.cellTypes,
        publication_citations: queryGroup2.publicationCitations,
      },
    };
  }, [organismId, queryGroups]);
}

export function useQueryGroupFilterDimensions(queryGroup: QueryGroup): {
  data: FilterDimensions;
  n_cells: number;
  isLoading: boolean;
} {
  const requestBody = useDEFiltersQueryRequestBodyForQueryGroups(queryGroup);

  const { data, isLoading } = useDEFiltersQuery(requestBody);

  return useMemo(() => {
    if (isLoading || !data)
      return { data: EMPTY_FILTER_DIMENSIONS, n_cells: 0, isLoading };

    const { filter_dims, n_cells } = data;

    const {
      datasets,
      development_stage_terms,
      disease_terms,
      self_reported_ethnicity_terms,
      sex_terms,
      tissue_terms,
      publication_citations,
      cell_type_terms,
      organism_terms,
    } = filter_dims;

    return {
      data: {
        datasets,
        development_stage_terms: development_stage_terms.map(toEntity),
        disease_terms: disease_terms.map(toEntity),
        self_reported_ethnicity_terms:
          self_reported_ethnicity_terms.map(toEntity),
        sex_terms: sex_terms.map(toEntity),
        publication_citations,
        tissue_terms: tissue_terms.map(toEntity),
        cell_type_terms: cell_type_terms.map(toEntity),
        organism_terms: organism_terms.map(toEntity),
      },
      n_cells,
      isLoading: false,
    };
  }, [data, isLoading]);
}

function useDEFiltersQueryRequestBodyForQueryGroups(queryGroup: QueryGroup) {
  const { organismId } = useContext(StateContext);
  const {
    developmentStages,
    diseases,
    ethnicities,
    sexes,
    tissues,
    cellTypes,
    publicationCitations,
  } = queryGroup;

  return useMemo(() => {
    if (!organismId) {
      return {
        filter: {
          development_stage_ontology_term_ids: [],
          disease_ontology_term_ids: [],
          organism_ontology_term_id: CXG_CENSUS_ORGANISMS[0].id,
          self_reported_ethnicity_ontology_term_ids: [],
          sex_ontology_term_ids: [],
          tissue_ontology_term_ids: [],
          cell_type_ontology_term_ids: [],
          publication_citations: [],
        },
      };
    }
    return {
      filter: {
        development_stage_ontology_term_ids: developmentStages,
        disease_ontology_term_ids: diseases,
        organism_ontology_term_id: organismId,
        self_reported_ethnicity_ontology_term_ids: ethnicities,
        sex_ontology_term_ids: sexes,
        tissue_ontology_term_ids: tissues,
        cell_type_ontology_term_ids: cellTypes,
        publication_citations: publicationCitations,
      },
    };
  }, [
    organismId,
    developmentStages,
    diseases,
    ethnicities,
    sexes,
    tissues,
    cellTypes,
    publicationCitations,
  ]);
}

function toEntity(item: RawOntologyTerm) {
  const [id, name] = Object.entries(item)[0];

  return { id, name: name || id || "" };
}

function useSnapshotId(): string | null {
  const state = useContext(StateContext);

  const { snapshotId } = state;

  return snapshotId || null;
}
