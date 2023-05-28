import { useMemo } from "react";
import { useQuery, UseQueryResult } from "react-query";
import { DEFAULT_FETCH_OPTIONS, JSON_BODY_FETCH_OPTIONS } from "./common";
import { ENTITIES } from "./entities";

export enum TYPES {
  CELL_ONTOLOGY_TREE = "CELL_ONTOLOGY_TREE",
  INITIAL_CELL_ONTOLOGY_TREE = "INITIAL_CELL_ONTOLOGY_TREE",
  SOURCE_DATA = "SOURCE_DATA",
  ENRICHED_GENES = "ENRICHED_GENES",
  CANONICAL_MARKERS = "CANONICAL_MARKERS",
  CL_DESCRIPTION = "CL_DESCRIPTION",
  DESCRIPTION = "DESCRIPTION",
  CELL_CARDS = "CELL_CARDS",
  INITIAL_CELL_ONTOLOGY_TREE_TISSUE = "INITIAL_CELL_ONTOLOGY_TREE_TISSUE",
  TISSUE_CARDS = "TISSUE_CARDS",
}

interface CellCardQuery {
  queryKey: {
    entities: ENTITIES[];
    id: string;
  };
  url: string;
}

export type CellCardResponse =
  | CellOntologyTreeResponse
  | InitialCellOntologyTreeStateResponse
  | SourceDataQueryResponse
  | EnrichedGenesQueryResponse
  | CanonicalMarkersQueryResponse
  | ClDescriptionQueryResponse
  | DescriptionQueryResponse
  | CellCardsQueryResponse;

/**
 * Generic fetch function
 */
async function fetchQuery({
  url,
  signal,
}: {
  url: string;
  signal: AbortSignal | undefined;
}): Promise<CellCardResponse | undefined> {
  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    method: "GET",
    signal,
  });
  if (response.status === 204) return undefined;
  const json: CellCardResponse = await response.json();

  if (!response.ok) {
    throw json;
  }

  return json;
}

/**
 * Generic cell cards hook
 */
export function useCellCardQuery<T = CellCardResponse>(
  dataType: TYPES,
  queryId = "" // Empty string if cell type is not needed for fetch function
): UseQueryResult<T> {
  const { queryKey, url: rawUrl } = QUERY_MAPPING[dataType];
  return useQuery(
    queryId ? [queryKey, queryId] : [queryKey],
    ({ signal }) =>
      fetchQuery({
        url: rawUrl.replace("%s", queryId), // Replacing raw url with entityId if applicable
        signal,
      }),
    {
      enabled: true,
      staleTime: Infinity,
    }
  );
}

/* ========== ontology_tree ========== */
export const USE_CELL_ONTOLOGY_TREE_QUERY = {
  entities: [ENTITIES.CELL_CARDS_CELL_ONTOLOGY_TREE],
  id: "cell-cards-cell-ontology-tree-query",
};

export interface CellOntologyTreeResponse {
  name: string;
  id: string;
  n_cells_rollup: number;
  n_cells: number;
  children?: this[];
}

export const useCellOntologyTree =
  (): UseQueryResult<CellOntologyTreeResponse> => {
    return useCellCardQuery<CellOntologyTreeResponse>(TYPES.CELL_ONTOLOGY_TREE);
  };

/* ========== ontology_tree_state ========== */
export const USE_INITIAL_CELL_ONTOLOGY_TREE_STATE_QUERY = {
  entities: [ENTITIES.CELL_CARDS_INITIAL_CELL_ONTOLOGY_TREE_STATE],
  id: "cell-cards-cell-ontology-tree-state-query",
};

export interface InitialCellOntologyTreeStateResponse {
  isExpandedNodes: string[];
  notShownWhenExpandedNodes: {
    [key: string]: string[];
  };
}

export const useCellOntologyTreeState = (
  entityId: string
): UseQueryResult<InitialCellOntologyTreeStateResponse> => {
  return useCellCardQuery<InitialCellOntologyTreeStateResponse>(
    TYPES.INITIAL_CELL_ONTOLOGY_TREE,
    entityId
  );
};

/* ========== ontology_tree_state_tissue ========== */
export const USE_INITIAL_CELL_ONTOLOGY_TREE_STATE_TISSUE_QUERY = {
  entities: [ENTITIES.CELL_CARDS_INITIAL_CELL_ONTOLOGY_TREE_STATE_TISSUE],
  id: "cell-cards-cell-ontology-tree-state-tissue-query",
};

export const useCellOntologyTreeStateTissue = (
  entityId: string
): UseQueryResult<InitialCellOntologyTreeStateResponse> => {
  return useCellCardQuery<InitialCellOntologyTreeStateResponse>(
    TYPES.INITIAL_CELL_ONTOLOGY_TREE_TISSUE,
    entityId
  );
};

/* ========== source_data ========== */
export const USE_SOURCE_DATA_QUERY = {
  entities: [ENTITIES.CELL_CARDS_SOURCE_DATA],
  id: "cell-cards-source-data-query",
};

interface SourceDataQueryResponseEntry {
  collection_name: string;
  collection_url: string;
  publication_url: string;
  publication_title: string;
  tissue: { label: string; ontology_term_id: string }[];
  disease: { label: string; ontology_term_id: string }[];
  organism: { label: string; ontology_term_id: string }[];
}

export type SourceDataQueryResponse = SourceDataQueryResponseEntry[];

export const useSourceData = (
  entityId: string
): UseQueryResult<SourceDataQueryResponse> => {
  return useCellCardQuery<SourceDataQueryResponse>(
    TYPES.SOURCE_DATA,
    entityId
  );
};

/* ========== enriched_genes ========== */
export const USE_ENRICHED_GENES_QUERY = {
  entities: [ENTITIES.CELL_CARDS_ENRICHED_GENES],
  id: "cell-cards-enriched-genes-query",
};

interface EnrichedGenesQueryResponseEntry {
  me: number;
  pc: number;
  symbol: string;
  name: string;
  organism: string;
}

export type EnrichedGenesQueryResponse = EnrichedGenesQueryResponseEntry[];

export const useEnrichedGenes = (
  entityId: string
): UseQueryResult<EnrichedGenesQueryResponse> => {
  return useCellCardQuery<EnrichedGenesQueryResponse>(
    TYPES.ENRICHED_GENES,
    entityId
  );
};

/* ========== canonical_markers ========== */
export const USE_CANONICAL_MARKERS_QUERY = {
  entities: [ENTITIES.CELL_CARDS_CANONICAL_MARKERS],
  id: "cell-cards-canonical-markersquery",
};

interface CanonicalMarkersQueryResponseEntry {
  tissue_general: string;
  tissue_specific: string;
  symbol: string;
  name: string;
  publication: string;
  publication_titles: string;
}

export type CanonicalMarkersQueryResponse =
  CanonicalMarkersQueryResponseEntry[];

export const useCanonicalMarkers = (
  entityId: string
): UseQueryResult<CanonicalMarkersQueryResponse> => {
  return useCellCardQuery<CanonicalMarkersQueryResponse>(
    TYPES.CANONICAL_MARKERS,
    entityId
  );
};

/* ========== CL description ========== */
export const USE_CL_DESCRIPTION_QUERY = {
  entities: [ENTITIES.CELL_CARDS_CL_DESCRIPTION],
  id: "cell-cards-cl-description-query",
};

export type ClDescriptionQueryResponse = string;

export const useClDescription = (
  entityId: string
): UseQueryResult<string> => {
  return useCellCardQuery<string>(TYPES.CL_DESCRIPTION, entityId);
};

/* ========== description ========== */
export const USE_DESCRIPTION_QUERY = {
  entities: [ENTITIES.CELL_CARDS_DESCRIPTION],
  id: "cell-cards-description-query",
};

export type DescriptionQueryResponse = string;

export const useDescription = (entityId: string): UseQueryResult<string> => {
  return useCellCardQuery<string>(TYPES.DESCRIPTION, entityId);
};

/* ========== cell_cards ========== */
export const USE_CELL_CARDS_QUERY = {
  entities: [ENTITIES.CELL_CARDS_CELL_CARDS],
  id: "cell-cards-query",
};

interface CellCardsQueryResponseEntry {
  id: string;
  label: string;
}

export type CellCardsQueryResponse = CellCardsQueryResponseEntry[];

export const useCellCards = (): UseQueryResult<CellCardsQueryResponse> => {
  return useCellCardQuery<CellCardsQueryResponse>(TYPES.CELL_CARDS);
};

/* ========== cell types by Id ========== */

export function useCellTypesById(): { [id: string]: string } | undefined {
  const { data, isLoading } = useCellCards();

  return useMemo(() => {
    if (!data || isLoading) return;
    const accumulator: { [id: string]: string } = {};
    return data.reduce((acc, curr) => {
      const { id, label } = curr;
      acc[id] = label;
      return acc;
    }, accumulator);
  }, [data, isLoading]);
}

/* ========== tissue_cards ========== */
export const USE_TISSUE_CARDS_QUERY = {
  entities: [ENTITIES.CELL_CARDS_TISSUE_CARDS],
  id: "tissue-cards-query",
};

interface TissueCardsQueryResponse {
  [tissueId: string]: {
    cell_types: string[];
    name: string;
  }
}

export const useTissueCards = (): UseQueryResult<TissueCardsQueryResponse> => {
  return useCellCardQuery<TissueCardsQueryResponse>(TYPES.TISSUE_CARDS);
};

/**
 * Mapping from data/response type to properties used for querying
 */
const QUERY_MAPPING: {
  [key in TYPES]: CellCardQuery;
} = {
  CELL_ONTOLOGY_TREE: {
    queryKey: USE_CELL_ONTOLOGY_TREE_QUERY,
    url: "/api/ontology_tree",
  },
  INITIAL_CELL_ONTOLOGY_TREE: {
    queryKey: USE_INITIAL_CELL_ONTOLOGY_TREE_STATE_QUERY,
    url: `/api/ontology_tree_state?entityId=%s`,
  },
  INITIAL_CELL_ONTOLOGY_TREE_TISSUE: {
    queryKey: USE_INITIAL_CELL_ONTOLOGY_TREE_STATE_TISSUE_QUERY,
    url: `/api/ontology_tree_state_tissue?entityId=%s`,
  },  
  SOURCE_DATA: {
    queryKey: USE_SOURCE_DATA_QUERY,
    url: `/api/source_data?entityId=%s`,
  },
  ENRICHED_GENES: {
    queryKey: USE_ENRICHED_GENES_QUERY,
    url: `/api/enriched_genes?entityId=%s`,
  },
  CANONICAL_MARKERS: {
    queryKey: USE_CANONICAL_MARKERS_QUERY,
    url: `/api/canonical_markers?entityId=%s`,
  },
  CL_DESCRIPTION: {
    queryKey: USE_CL_DESCRIPTION_QUERY,
    url: `/api/cl_description?entityId=%s`,
  },
  DESCRIPTION: {
    queryKey: USE_DESCRIPTION_QUERY,
    url: `/api/description?entityId=%s`,
  },
  CELL_CARDS: {
    queryKey: USE_CELL_CARDS_QUERY,
    url: "/api/cell_cards",
  },
  TISSUE_CARDS: {
    queryKey: USE_TISSUE_CARDS_QUERY,
    url: "/api/tissue_cards",
  },  
};
