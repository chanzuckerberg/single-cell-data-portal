import { useContext, useMemo } from "react";
import { useQuery, UseQueryResult } from "react-query";
import { API_URL } from "src/configs/dev";
import {
  DispatchContext,
  StateContext,
} from "src/views/DifferentialExpression/common/store";
import { setSnapshotId } from "src/views/DifferentialExpression/common/store/actions";
import { API } from "../API";
import { ROUTES } from "../constants/routes";
import { DEFAULT_FETCH_OPTIONS, JSON_BODY_FETCH_OPTIONS } from "./common";
import { ENTITIES } from "./entities";
import {
  EMPTY_FILTERS,
  QueryGroup,
} from "src/views/DifferentialExpression/common/store/reducer";

interface RawOntologyTerm {
  [id: string]: string;
}

interface RawOntologyTermsByOrganism {
  [organismId: string]: Array<RawOntologyTerm>;
}
export interface RawPrimaryFilterDimensionsResponse {
  gene_terms: RawOntologyTermsByOrganism;
  organism_terms: Array<RawOntologyTerm>;
  snapshot_id: string;
  tissue_terms: RawOntologyTermsByOrganism;
}

interface FilterSecondary {
  organism_ontology_term_id: string;
  tissue_ontology_term_ids: string[];
  dataset_ids: string[];
  disease_ontology_term_ids: string[];
  sex_ontology_term_ids: string[];
  development_stage_ontology_term_ids: string[];
  self_reported_ethnicity_ontology_term_ids: string[];
  cell_type_ontology_term_ids: string[];
}

export interface FiltersQuery {
  filter: FilterSecondary;
}

export interface DifferentialExpressionQuery {
  queryGroup1Filters: {
    dataset_ids: string[];
    development_stage_ontology_term_ids: string[];
    disease_ontology_term_ids: string[];
    organism_ontology_term_id: string;
    self_reported_ethnicity_ontology_term_ids: string[];
    sex_ontology_term_ids: string[];
    tissue_ontology_term_ids: string[];
    cell_type_ontology_term_ids: string[];
  };
  queryGroup2Filters: {
    dataset_ids: string[];
    development_stage_ontology_term_ids: string[];
    disease_ontology_term_ids: string[];
    organism_ontology_term_id: string;
    self_reported_ethnicity_ontology_term_ids: string[];
    sex_ontology_term_ids: string[];
    tissue_ontology_term_ids: string[];
    cell_type_ontology_term_ids: string[];
  };
}

interface FiltersQueryResponse {
  filter_dims: {
    datasets: {
      collection_id: string;
      collection_label: string;
      id: string;
      label: string;
    }[];
    disease_terms: { [id: string]: string }[];
    sex_terms: { [id: string]: string }[];
    development_stage_terms: { [id: string]: string }[];
    self_reported_ethnicity_terms: { [id: string]: string }[];
    tissue_terms: { [id: string]: string }[];
    cell_type_terms: { [id: string]: string }[];
    organism_terms: { [id: string]: string }[];
  };
  snapshot_id: string;
}

interface GetDeQueryResponse {
  query_criteria1: {
    disease_terms: { [id: string]: string }[];
    sex_terms: { [id: string]: string }[];
    development_stage_terms: { [id: string]: string }[];
    self_reported_ethnicity_terms: { [id: string]: string }[];
    tissue_terms: { [id: string]: string }[];
    cell_type_terms: { [id: string]: string }[];
    organism_terms: { [id: string]: string }[];
  };
  query_criteria2: {
    disease_terms: { [id: string]: string }[];
    sex_terms: { [id: string]: string }[];
    development_stage_terms: { [id: string]: string }[];
    self_reported_ethnicity_terms: { [id: string]: string }[];
    tissue_terms: { [id: string]: string }[];
    cell_type_terms: { [id: string]: string }[];
    organism_terms: { [id: string]: string }[];
  };
  snapshot_id: string;
}

export interface DifferentialExpressionResult {
  gene_ontology_term_id: string;
  gene_symbol: string;
  p_value: number;
  effect_size: number;
}

export interface PathwayEnrichmentAnalysisResult {
  gene_set: string;
  gene_symbols: string;
  p_value: number;
  fdr_q_value: number;
}

interface DifferentialExpressionQueryResponse {
  differentialExpressionResults1: DifferentialExpressionResult[];
  differentialExpressionResults2: DifferentialExpressionResult[];
  pathwayEnrichmentAnalysisResults1: PathwayEnrichmentAnalysisResult[];
  pathwayEnrichmentAnalysisResults2: PathwayEnrichmentAnalysisResult[];
  successCode: number;
  snapshot_id: string;
}

interface DifferentialExpressionQueryResult {
  differentialExpressionResults1: DifferentialExpressionResult[];
  differentialExpressionResults2: DifferentialExpressionResult[];
  pathwayEnrichmentAnalysisResults1: PathwayEnrichmentAnalysisResult[];
  pathwayEnrichmentAnalysisResults2: PathwayEnrichmentAnalysisResult[];
  successCode?: number;
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

export const USE_FILTERS_QUERY = {
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

async function fetchGetDeQuery({
  query,
  signal,
}: {
  query: { user_query: string };
  signal?: AbortSignal;
}): Promise<GetDeQueryResponse | undefined> {
  if (!query) return;

  const url = API_URL + API.DE_GET_QUERY;

  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    body: JSON.stringify(query),
    method: "POST",
    signal,
  });
  const json: GetDeQueryResponse = await response.json();
  if (!response.ok) {
    throw json;
  }

  return json;
}

export const USE_DE_QUERY = {
  entities: [ENTITIES.DE_QUERY],
  id: "de-query",
};

export const USE_GET_DE_QUERY = {
  entities: [ENTITIES.GET_DE_QUERY],
  id: "get-de-query",
};

export function useWMGFiltersQuery(
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

export function useDEQuery(
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
  development_stage_terms: [],
  disease_terms: [],
  self_reported_ethnicity_terms: [],
  sex_terms: [],
  tissue_terms: [],
  cell_type_terms: [],
  organism_terms: [],
};

export interface RawDataset {
  collection_id: string;
  collection_label: string;
  id: string;
  label: string;
}

export interface FilterDimensions {
  datasets: RawDataset[];
  development_stage_terms: { id: string; name: string }[];
  disease_terms: { id: string; name: string }[];
  self_reported_ethnicity_terms: { id: string; name: string }[];
  sex_terms: { id: string; name: string }[];
  tissue_terms: { id: string; name: string }[];
  cell_type_terms: { id: string; name: string }[];
  organism_terms: { id: string; name: string }[];
}

const TEMP_ALLOW_NAME_LIST = ["Homo sapiens", "Mus musculus"];

export function useAvailableOrganisms() {
  const { data, isLoading } = useQueryGroupFilterDimensions(EMPTY_FILTERS);

  return useMemo(() => {
    if (isLoading || !data) return { data: [], isLoading };

    const { organism_terms } = data;

    return {
      isLoading,
      data: organism_terms.filter((organism) =>
        TEMP_ALLOW_NAME_LIST.includes(organism.name)
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
    if (isLoading || !data || (data?.successCode ?? 0) > 0)
      return {
        data: {
          differentialExpressionResults1: [],
          differentialExpressionResults2: [],
          pathwayEnrichmentAnalysisResults1: [],
          pathwayEnrichmentAnalysisResults2: [],
          successCode: data?.successCode ?? 0,
        },
        isLoading,
      };
    return {
      data: {
        differentialExpressionResults1: data.differentialExpressionResults1,
        differentialExpressionResults2: data.differentialExpressionResults2,
        pathwayEnrichmentAnalysisResults1:
          data.pathwayEnrichmentAnalysisResults1,
        pathwayEnrichmentAnalysisResults2:
          data.pathwayEnrichmentAnalysisResults2,
        successCode: 0,
      },
      isLoading: false,
    };
  }, [data, isLoading]);
}

function useGetDeQueryRequestBody(userInput: string) {
  return {
    user_query: userInput,
  };
}

export function useGetDeQuery(query: {
  user_query: string;
}): UseQueryResult<GetDeQueryResponse> {
  const dispatch = useContext(DispatchContext);

  // (thuang): Refresh query when the snapshotId changes
  const currentSnapshotId = useSnapshotId();

  return useQuery(
    [USE_GET_DE_QUERY, query, currentSnapshotId],
    ({ signal }) => fetchGetDeQuery({ query, signal }),
    {
      enabled: query.user_query !== "",
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

interface NaturalLanguageDeQuery {
  organism: string;
  queryCriteria1: QueryGroup;
  queryCriteria2: QueryGroup;
  queryCriteriaNames1: QueryGroup;
  queryCriteriaNames2: QueryGroup;
  isLoading: boolean;
}

export function useNaturalLanguageDeQuery(
  userInput: string
): NaturalLanguageDeQuery {
  const { organismId } = useContext(StateContext);
  const requestBody = useGetDeQueryRequestBody(userInput);
  const { data, isLoading } = useGetDeQuery(requestBody);
  return useMemo(() => {
    if (isLoading || !data || !organismId)
      return {
        organism: "",
        queryCriteria1: EMPTY_FILTERS,
        queryCriteria2: EMPTY_FILTERS,
        queryCriteriaNames1: EMPTY_FILTERS,
        queryCriteriaNames2: EMPTY_FILTERS,
        isLoading,
      };
    let organism = organismId;
    if (data.query_criteria1.organism_terms) {
      organism = data.query_criteria1.organism_terms.map(toId)[0];
    } else if (data.query_criteria2.organism_terms) {
      organism = data.query_criteria2.organism_terms.map(toId)[0];
    }
    const queryCriteria1 = _formatQueryCriteria(data.query_criteria1, toId);
    const queryCriteria2 = _formatQueryCriteria(data.query_criteria2, toId);
    const queryCriteriaNames1 = _formatQueryCriteria(
      data.query_criteria1,
      toName
    );
    const queryCriteriaNames2 = _formatQueryCriteria(
      data.query_criteria2,
      toName
    );

    return {
      organism,
      queryCriteria1,
      queryCriteria2,
      queryCriteriaNames1,
      queryCriteriaNames2,
      isLoading: false,
    };
  }, [data, isLoading, organismId]);
}

function _formatQueryCriteria(
  queryCriteria: GetDeQueryResponse["query_criteria1"],
  mapper: (item: RawOntologyTerm) => string
) {
  return {
    datasets: [],
    diseases: queryCriteria.disease_terms?.map(mapper) ?? [],
    cellTypes: queryCriteria.cell_type_terms?.map(mapper) ?? [],
    developmentStages: queryCriteria.development_stage_terms?.map(mapper) ?? [],
    ethnicities: queryCriteria.self_reported_ethnicity_terms?.map(mapper) ?? [],
    sexes: queryCriteria.sex_terms?.map(mapper) ?? [],
    tissues: queryCriteria.tissue_terms?.map(mapper) ?? [],
  };
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
        dataset_ids: queryGroup1.datasets,
        development_stage_ontology_term_ids: queryGroup1.developmentStages,
        disease_ontology_term_ids: queryGroup1.diseases,
        organism_ontology_term_id: organismId,
        self_reported_ethnicity_ontology_term_ids: queryGroup1.ethnicities,
        sex_ontology_term_ids: queryGroup1.sexes,
        tissue_ontology_term_ids: queryGroup1.tissues,
        cell_type_ontology_term_ids: queryGroup1.cellTypes,
      },
      queryGroup2Filters: {
        dataset_ids: queryGroup2.datasets,
        development_stage_ontology_term_ids: queryGroup2.developmentStages,
        disease_ontology_term_ids: queryGroup2.diseases,
        organism_ontology_term_id: organismId,
        self_reported_ethnicity_ontology_term_ids: queryGroup2.ethnicities,
        sex_ontology_term_ids: queryGroup2.sexes,
        tissue_ontology_term_ids: queryGroup2.tissues,
        cell_type_ontology_term_ids: queryGroup2.cellTypes,
      },
    };
  }, [organismId, queryGroups]);
}

export function useQueryGroupFilterDimensions(queryGroup: QueryGroup): {
  data: FilterDimensions;
  isLoading: boolean;
} {
  const requestBody = useWMGFiltersQueryRequestBodyForQueryGroups(queryGroup);
  const { data, isLoading } = useWMGFiltersQuery(requestBody);

  return useMemo(() => {
    if (isLoading || !data) return { data: EMPTY_FILTER_DIMENSIONS, isLoading };

    const { filter_dims } = data;

    const {
      datasets,
      development_stage_terms,
      disease_terms,
      self_reported_ethnicity_terms,
      sex_terms,
      tissue_terms,
      cell_type_terms,
      organism_terms,
    } = filter_dims;

    const sortedDatasets = Object.values(
      aggregateCollectionsFromDatasets(datasets)
    ).flatMap(({ datasets }) => datasets);
    return {
      data: {
        datasets: sortedDatasets.map((dataset) => ({
          ...dataset,
          name: dataset.label,
        })),
        development_stage_terms: development_stage_terms.map(toEntity),
        disease_terms: disease_terms.map(toEntity),
        self_reported_ethnicity_terms:
          self_reported_ethnicity_terms.map(toEntity),
        sex_terms: sex_terms.map(toEntity),
        tissue_terms: tissue_terms.map(toEntity),
        cell_type_terms: cell_type_terms.map(toEntity),
        organism_terms: organism_terms.map(toEntity),
      },
      isLoading: false,
    };
  }, [data, isLoading]);
}

function useWMGFiltersQueryRequestBodyForQueryGroups(queryGroup: QueryGroup) {
  const { organismId } = useContext(StateContext);
  const {
    datasets,
    developmentStages,
    diseases,
    ethnicities,
    sexes,
    tissues,
    cellTypes,
  } = queryGroup;

  return useMemo(() => {
    if (!organismId) {
      return {
        filter: {
          dataset_ids: [],
          development_stage_ontology_term_ids: [],
          disease_ontology_term_ids: [],
          organism_ontology_term_id: TEMP_ALLOW_NAME_LIST[0],
          self_reported_ethnicity_ontology_term_ids: [],
          sex_ontology_term_ids: [],
          tissue_ontology_term_ids: [],
          cell_type_ontology_term_ids: [],
        },
      };
    }
    return {
      filter: {
        dataset_ids: datasets,
        development_stage_ontology_term_ids: developmentStages,
        disease_ontology_term_ids: diseases,
        organism_ontology_term_id: organismId,
        self_reported_ethnicity_ontology_term_ids: ethnicities,
        sex_ontology_term_ids: sexes,
        tissue_ontology_term_ids: tissues,
        cell_type_ontology_term_ids: cellTypes,
      },
    };
  }, [
    organismId,
    datasets,
    developmentStages,
    diseases,
    ethnicities,
    sexes,
    tissues,
    cellTypes,
  ]);
}

function toEntity(item: RawOntologyTerm) {
  const [id, name] = Object.entries(item)[0];

  return { id, name: name || id || "" };
}

function toId(item: RawOntologyTerm) {
  return Object.keys(item)[0];
}
function toName(item: RawOntologyTerm) {
  return Object.values(item)[0];
}

function useSnapshotId(): string | null {
  const state = useContext(StateContext);

  const { snapshotId } = state;

  return snapshotId || null;
}

interface Dataset extends RawDataset {
  id: string;
  label: string;
}

export interface CollectionFromDatasets {
  name: string;
  url: string;
  datasets: Dataset[];
}

export interface CollectionsFromDatasets {
  [name: string]: CollectionFromDatasets;
}

export function aggregateCollectionsFromDatasets(
  datasets: FilterDimensions["datasets"]
): CollectionsFromDatasets {
  const collections: CollectionsFromDatasets = {};

  for (const dataset of datasets) {
    const { collection_label, collection_id, id, label } = dataset;

    if (!collections[collection_label]) {
      collections[collection_label] = {
        datasets: [],
        name: collection_label,
        url: ROUTES.COLLECTION.replace(":id", collection_id),
      };
    }

    collections[collection_label].datasets.push({ ...dataset, id, label });
  }

  for (const collection of Object.values(collections)) {
    collection.datasets.sort((a, b) => {
      return a.label.localeCompare(b.label);
    });
  }

  return collections;
}
