import { useContext, useMemo } from "react";
import { useQuery, UseQueryResult } from "react-query";
import { API_URL } from "src/configs/configs";
import { State, StateContext } from "src/views/WheresMyGene/common/store";
import {
  CellType,
  CellTypeGeneExpressionSummaryData,
  GeneExpressionSummary,
  RawCellTypeGeneExpressionSummaryData,
  Tissue,
} from "src/views/WheresMyGene/common/types";
import { API } from "../API";
import { EMPTY_OBJECT } from "../constants/utils";
import { DEFAULT_FETCH_OPTIONS, JSON_BODY_FETCH_OPTIONS } from "./common";
import { CELL_TYPE_ORDER_BY_TISSUE } from "./constants/cellTypeOrderByTissue";
import { ENTITIES } from "./entities";

interface RawPrimaryFilterDimensionsResponse {
  gene_terms: { [organismId: string]: Array<{ [id: string]: string }> };
  organism_terms: { [id: string]: string }[];
  snapshot_id: string;
  tissue_terms: { [id: string]: string }[];
}

export interface PrimaryFilterDimensionsResponse {
  genes: {
    [organismId: string]: {
      id: string;
      name: string;
    }[];
  };
  organisms: {
    id: string;
    name: string;
  }[];
  snapshotId: string;
  tissues: {
    id: string;
    name: string;
  }[];
}

export async function fetchPrimaryFilterDimensions(): Promise<PrimaryFilterDimensionsResponse> {
  const url = API_URL + API.WMG_PRIMARY_FILTER_DIMENSIONS;

  const response = await (await fetch(url, DEFAULT_FETCH_OPTIONS)).json();

  return transformPrimaryFilterDimensions(response);
}

function transformPrimaryFilterDimensions(
  response: RawPrimaryFilterDimensionsResponse
): PrimaryFilterDimensionsResponse {
  const { gene_terms, organism_terms, snapshot_id, tissue_terms } = response;

  return {
    genes: Object.entries(gene_terms).reduce((memo, [organismId, genes]) => {
      memo[organismId] = genes.map(toEntity);
      return memo;
    }, {} as { [organismId: string]: { id: string; name: string }[] }),
    organisms: organism_terms.map(toEntity),
    snapshotId: snapshot_id,
    tissues: tissue_terms.map(toEntity),
  };
}

export const USE_PRIMARY_FILTER_DIMENSIONS = {
  entities: [ENTITIES.WMG_PRIMARY_FILTER_DIMENSIONS],
  id: "wmg-primaryFilterDimensions",
};

export function usePrimaryFilterDimensions(): UseQueryResult<PrimaryFilterDimensionsResponse> {
  return useQuery<PrimaryFilterDimensionsResponse>(
    [USE_PRIMARY_FILTER_DIMENSIONS],
    fetchPrimaryFilterDimensions,
    // (thuang): We don't need to refetch during the session
    { staleTime: Infinity }
  );
}

interface Filter {
  gene_ontology_term_ids: string[];
  organism_ontology_term_id: string;
  tissue_ontology_term_ids: string[];
  dataset_ids?: string[];
  disease_ontology_term_ids?: string[];
  sex_ontology_term_ids?: string[];
  development_stage_ontology_term_ids?: string[];
  ethnicity_ontology_term_ids?: string[];
}

export interface Query {
  snapshot_id: string;
  include_filter_dims: boolean;
  filter: Filter;
}

interface QueryResponse {
  expression_summary: {
    // gene_ontology_term_id
    [geneId: string]: {
      [tissueId: string]: {
        // cell_type_ontology_term_id
        id: string;
        me: number;
        n: number;
        pc: number;
        tpc: number;
      }[];
    };
  };
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
    ethnicity_terms: { [id: string]: string }[];
  };
  snapshot_id: string;
  term_id_labels: {
    cell_types: {
      [tissue_type_ontology_term_id: string]: {
        [id: string]: string;
      }[];
    };
    genes: {
      [id: string]: string;
    }[];
  };
}

async function fetchQuery(
  query: Query | null
): Promise<QueryResponse | undefined> {
  if (!query) return;

  const url = API_URL + API.WMG_QUERY;
  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    body: JSON.stringify(query),
    method: "POST",
  });
  const json = await response.json();

  if (!response.ok) {
    throw json;
  }

  return json;
}

export const USE_QUERY = {
  entities: [ENTITIES.WMG_QUERY],
  id: "wmg-query",
};

export function useWMGQuery(
  query: Query | null
): UseQueryResult<QueryResponse> {
  return useQuery([USE_QUERY, query], () => fetchQuery(query), {
    enabled: Boolean(query),
    // (thuang): We don't need to refetch during the session
    staleTime: Infinity,
  });
}

const EMPTY_FILTER_DIMENSIONS = {
  datasets: [],
  development_stage_terms: [],
  disease_terms: [],
  ethnicity_terms: [],
  sex_terms: [],
};

export interface FilterDimensions {
  datasets: {
    collection_id: string;
    collection_label: string;
    id: string;
    label: string;
  }[];
  development_stage_terms: { id: string; name: string }[];
  disease_terms: { id: string; name: string }[];
  ethnicity_terms: { id: string; name: string }[];
  sex_terms: { id: string; name: string }[];
}

/**
 * (thuang): For Filters panel, `includeAllFilterOptions` should be `true`, so BE
 * returns all available secondary filter options for us to display
 */
export function useFilterDimensions(
  options = { includeAllFilterOptions: false }
): {
  data: FilterDimensions;
  isLoading: boolean;
} {
  const { includeAllFilterOptions } = options;

  const requestBody = useWMGQueryRequestBody({ includeAllFilterOptions });

  const { data, isLoading } = useWMGQuery(requestBody);

  return useMemo(() => {
    if (isLoading || !data) return { data: EMPTY_FILTER_DIMENSIONS, isLoading };

    const { filter_dims } = data;

    const {
      datasets,
      development_stage_terms,
      disease_terms,
      ethnicity_terms,
      sex_terms,
    } = filter_dims;

    return {
      data: {
        datasets: datasets.map((dataset) => ({
          ...dataset,
          name: dataset.label,
        })),
        development_stage_terms: development_stage_terms.map(toEntity),
        disease_terms: disease_terms.map(toEntity),
        ethnicity_terms: ethnicity_terms.map(toEntity),
        sex_terms: sex_terms.map(toEntity),
      },
      isLoading: false,
    };
  }, [data, isLoading]);
}

export function useExpressionSummary(): {
  isLoading: boolean;
  data: QueryResponse["expression_summary"];
} {
  const requestBody = useWMGQueryRequestBody();

  const { data, isLoading } = useWMGQuery(requestBody);

  return useMemo(() => {
    if (isLoading || !data) return { data: EMPTY_OBJECT, isLoading };

    const { expression_summary } = data;

    return {
      data: expression_summary,
      isLoading: false,
    };
  }, [data, isLoading]);
}

export interface CellTypeByTissueName {
  [tissueName: string]: CellType[];
}

export function useCellTypesByTissueName(): {
  isLoading: boolean;
  data: CellTypeByTissueName;
} {
  const { data, isLoading } = useExpressionSummary();

  const {
    data: primaryFilterDimensions,
    isLoading: isLoadingPrimaryFilterDimensions,
  } = usePrimaryFilterDimensions();

  const { data: termIdLabels, isLoading: isLoadingTermIdLabels } =
    useTermIdLabels();

  return useMemo(() => {
    if (
      isLoading ||
      !data ||
      isLoadingPrimaryFilterDimensions ||
      !primaryFilterDimensions ||
      isLoadingTermIdLabels ||
      !termIdLabels
    ) {
      return { data: EMPTY_OBJECT, isLoading };
    }

    const result: Map<Tissue, Map<string, CellType>> = new Map();
    const { tissues } = primaryFilterDimensions;

    const tissuesById: { [id: string]: { id: string; name: string } } = {};

    for (const tissue of tissues) {
      tissuesById[tissue.id] = tissue;
    }

    for (const expressionSummaryByTissue of Object.values(data)) {
      for (const [tissueID, expressionSummaries] of Object.entries(
        expressionSummaryByTissue
      )) {
        const cellTypes = result.get(tissueID) || new Map();

        for (const expressionSummary of expressionSummaries) {
          const { id } = expressionSummary;

          const cellType = {
            id,
            name: termIdLabels.cell_types[tissueID][id],
          };

          cellTypes.set(id, cellType);
        }

        result.set(tissueID, cellTypes);
      }
    }

    const cellTypesByTissueName: { [tissueName: string]: CellType[] } = {};

    for (const [tissueId, cellTypesById] of result.entries()) {
      const tissueName = tissuesById[tissueId].name;

      const cellTypeOrder = CELL_TYPE_ORDER_BY_TISSUE[tissueId];

      const cellTypes = Array.from(cellTypesById.values());

      if (cellTypeOrder) {
        cellTypesByTissueName[tissueName] = cellTypes.sort((a, b) => {
          const aIndex = cellTypeOrder[a.id];
          const bIndex = cellTypeOrder[b.id];

          return aIndex - bIndex;
        });
      } else {
        console.warn("No cell type order for tissue", tissueId);
        cellTypesByTissueName[tissueName] = cellTypes;
      }
    }

    return {
      data: cellTypesByTissueName,
      isLoading,
    };
  }, [
    data,
    isLoading,
    primaryFilterDimensions,
    isLoadingPrimaryFilterDimensions,
    termIdLabels,
    isLoadingTermIdLabels,
  ]);
}

export interface GeneExpressionSummariesByTissueName {
  [tissueName: string]: { [geneName: string]: GeneExpressionSummary };
}

export function useGeneExpressionSummariesByTissueName(): {
  data: GeneExpressionSummariesByTissueName;
  isLoading: boolean;
} {
  const { data, isLoading } = useExpressionSummary();

  const {
    data: primaryFilterDimensions,
    isLoading: isLoadingPrimaryFilterDimensions,
  } = usePrimaryFilterDimensions();

  const { data: termIdLabels, isLoading: isLoadingTermIdLabels } =
    useTermIdLabels();

  return useMemo(() => {
    if (
      isLoading ||
      !data ||
      isLoadingPrimaryFilterDimensions ||
      !primaryFilterDimensions ||
      isLoadingTermIdLabels ||
      !termIdLabels
    ) {
      return { data: EMPTY_OBJECT, isLoading };
    }

    const { tissues } = primaryFilterDimensions;

    const tissuesById: { [id: string]: { id: string; name: string } } = {};

    for (const tissue of tissues) {
      tissuesById[tissue.id] = tissue;
    }

    const result: {
      [tissueName: string]: { [geneName: string]: GeneExpressionSummary };
    } = {};

    for (const [geneId, expressionSummariesByTissue] of Object.entries(data)) {
      const geneName = termIdLabels.genes[geneId];

      for (const [tissueId, expressionSummaries] of Object.entries(
        expressionSummariesByTissue
      )) {
        const tissueName = tissuesById[tissueId].name;

        const tissueGeneExpressionSummaries: {
          [geneName: string]: GeneExpressionSummary;
        } = result[tissueName] || {};

        tissueGeneExpressionSummaries[geneName] = {
          cellTypeGeneExpressionSummaries: expressionSummaries.map(
            transformCellTypeGeneExpressionSummaryData
          ),
          name: geneName,
        };

        result[tissueName] = tissueGeneExpressionSummaries;
      }
    }

    return { data: result, isLoading };
  }, [
    data,
    isLoading,
    primaryFilterDimensions,
    isLoadingPrimaryFilterDimensions,
    termIdLabels,
    isLoadingTermIdLabels,
  ]);
}

function transformCellTypeGeneExpressionSummaryData(
  data: RawCellTypeGeneExpressionSummaryData
): CellTypeGeneExpressionSummaryData {
  const { id, pc, me } = data;

  return {
    ...data,
    id,
    meanExpression: me,
    percentage: pc,
  };
}

interface TermIdLabels {
  cell_types: { [tissueID: string]: { [id: string]: string } };
  genes: { [id: string]: string };
}

export function useTermIdLabels(): {
  data: TermIdLabels;
  isLoading: boolean;
} {
  const requestBody = useWMGQueryRequestBody();
  const { data, isLoading } = useWMGQuery(requestBody);

  return useMemo(() => {
    if (isLoading || !data)
      return {
        data: { cell_types: EMPTY_OBJECT, genes: EMPTY_OBJECT },
        isLoading,
      };

    const {
      term_id_labels: { cell_types, genes },
    } = data;

    const returnCellTypes: TermIdLabels["cell_types"] = {};
    Object.entries(cell_types).forEach(([tissueID, cell_types]) => {
      returnCellTypes[tissueID] = aggregateIdLabels(cell_types);
    });

    return {
      data: {
        cell_types: returnCellTypes,
        genes: aggregateIdLabels(genes),
      },
      isLoading: false,
    };
  }, [data, isLoading]);
}

function aggregateIdLabels(items: { [id: string]: string }[]): {
  [id: string]: string;
} {
  return items.reduce((memo, item) => ({ ...memo, ...item }), {});
}

const EMPTY_FILTERS: State["selectedFilters"] = {
  datasets: undefined,
  developmentStages: undefined,
  diseases: undefined,
  ethnicities: undefined,
  sexes: undefined,
};

function useWMGQueryRequestBody(options = { includeAllFilterOptions: false }) {
  const { includeAllFilterOptions } = options;

  const {
    selectedGenes,
    selectedTissues,
    selectedOrganismId,
    selectedFilters,
  } = useContext(StateContext);

  const { data } = usePrimaryFilterDimensions();

  /**
   * (thuang): When `includeAllFilterOptions` is `true`, we don't want to pass
   * any selected secondary filter options to the query, otherwise BE will return
   * only the filtered options back to us.
   */
  const { datasets, developmentStages, diseases, ethnicities, sexes } =
    includeAllFilterOptions ? EMPTY_FILTERS : selectedFilters;

  const organismGenesByName = useMemo(() => {
    const result: { [name: string]: { id: string; name: string } } = {};

    if (!data || !selectedOrganismId) return result;

    const { genes } = data;

    const organismGenes = genes[selectedOrganismId];

    for (const gene of organismGenes) {
      result[gene.name] = gene;
    }

    return result;
  }, [data, selectedOrganismId]);

  const tissuesByName = useMemo(() => {
    const result: { [name: string]: { id: string; name: string } } = {};

    if (!data) return result;

    const { tissues } = data;

    for (const tissue of tissues) {
      result[tissue.name] = tissue;
    }

    return result;
  }, [data]);

  return useMemo(() => {
    if (
      !data ||
      !selectedOrganismId ||
      !selectedTissues.length ||
      !selectedGenes.length
    ) {
      return null;
    }

    const { snapshotId } = data;

    const gene_ontology_term_ids = selectedGenes.map((geneName) => {
      return organismGenesByName[geneName].id;
    });

    const tissue_ontology_term_ids = selectedTissues.map((tissueName) => {
      return tissuesByName[tissueName].id;
    });

    return {
      filter: {
        dataset_ids: datasets,
        development_stage_ontology_term_ids: developmentStages,
        disease_ontology_term_ids: diseases,
        ethnicity_ontology_term_ids: ethnicities,
        gene_ontology_term_ids,
        organism_ontology_term_id: selectedOrganismId,
        sex_ontology_term_ids: sexes,
        tissue_ontology_term_ids,
      },
      include_filter_dims: true,
      snapshot_id: snapshotId,
    };
  }, [
    selectedGenes,
    selectedTissues,
    selectedOrganismId,
    data,
    organismGenesByName,
    tissuesByName,
    datasets,
    developmentStages,
    diseases,
    ethnicities,
    sexes,
  ]);
}

function toEntity(item: { [id: string]: string }) {
  const [id, name] = Object.entries(item)[0];

  return { id, name };
}
