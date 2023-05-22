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
  cellTypeId = "" // Empty string if cell type is not needed for fetch function
): UseQueryResult<T> {
  const { queryKey, url: rawUrl } = QUERY_MAPPING[dataType];

  return useQuery(
    cellTypeId ? [queryKey, cellTypeId] : [queryKey],
    ({ signal }) =>
      fetchQuery({
        url: rawUrl.replace("%s", cellTypeId), // Replacing raw url with cellTypeId if applicable
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
  id: "cell-explorer-cell-ontology-tree-query",
};

export interface CellOntologyTreeResponse {
  name: string;
  id: string;
  n_cells_rollup: number;
  n_cells: number;
  children?: this[];
  _children?: this[];
}

export const useCellOntologyTree =
  (): UseQueryResult<CellOntologyTreeResponse> => {
    return useCellCardQuery<CellOntologyTreeResponse>(TYPES.CELL_ONTOLOGY_TREE);
  };

/* ========== ontology_tree_state ========== */
export const USE_INITIAL_CELL_ONTOLOGY_TREE_STATE_QUERY = {
  entities: [ENTITIES.CELL_CARDS_INITIAL_CELL_ONTOLOGY_TREE_STATE],
  id: "cell-explorer-cell-ontology-tree-query",
};

export interface InitialCellOntologyTreeStateResponse {
  isExpandedNodes: string[];
  notShownWhenExpandedNodes: {
    [key: string]: string[];
  };
}

export const useCellOntologyTreeState = (
  cellTypeId: string
): UseQueryResult<InitialCellOntologyTreeStateResponse> => {
  return useCellCardQuery<InitialCellOntologyTreeStateResponse>(
    TYPES.INITIAL_CELL_ONTOLOGY_TREE,
    cellTypeId
  );
};

/* ========== source_data ========== */
export const USE_SOURCE_DATA_QUERY = {
  entities: [ENTITIES.CELL_CARDS_SOURCE_DATA],
  id: "cell-explorer-source-data-query",
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
  cellTypeId: string
): UseQueryResult<SourceDataQueryResponse> => {
  return useCellCardQuery<SourceDataQueryResponse>(
    TYPES.SOURCE_DATA,
    cellTypeId
  );
};

/* ========== enriched_genes ========== */
export const USE_ENRICHED_GENES_QUERY = {
  entities: [ENTITIES.CELL_CARDS_ENRICHED_GENES],
  id: "cell-explorer-enriched-genes-query",
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
  cellTypeId: string
): UseQueryResult<EnrichedGenesQueryResponse> => {
  return useCellCardQuery<EnrichedGenesQueryResponse>(
    TYPES.ENRICHED_GENES,
    cellTypeId
  );
};

/* ========== canonical_markers ========== */
export const USE_CANONICAL_MARKERS_QUERY = {
  entities: [ENTITIES.CELL_CARDS_CANONICAL_MARKERS],
  id: "cell-explorer-canonical-markersquery",
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
  cellTypeId: string
): UseQueryResult<CanonicalMarkersQueryResponse> => {
  return useCellCardQuery<CanonicalMarkersQueryResponse>(
    TYPES.CANONICAL_MARKERS,
    cellTypeId
  );
};

/* ========== CL description ========== */
export const USE_CL_DESCRIPTION_QUERY = {
  entities: [ENTITIES.CELL_CARDS_CL_DESCRIPTION],
  id: "cell-explorer-cl-description-query",
};

export type ClDescriptionQueryResponse = string;

export const useClDescription = (
  cellTypeId: string
): UseQueryResult<string> => {
  return useCellCardQuery<string>(TYPES.CL_DESCRIPTION, cellTypeId);
};

/* ========== description ========== */
export const USE_DESCRIPTION_QUERY = {
  entities: [ENTITIES.CELL_CARDS_DESCRIPTION],
  id: "cell-explorer-description-query",
};

export type DescriptionQueryResponse = string;

export const useDescription = (cellTypeId: string): UseQueryResult<string> => {
  return useCellCardQuery<string>(TYPES.DESCRIPTION, cellTypeId);
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

export const useCellTypes = (): UseQueryResult<CellCardsQueryResponse> => {
  return useCellCardQuery<CellCardsQueryResponse>(TYPES.CELL_CARDS);
};

/* ========== cell types by Id ========== */

export function useCellTypesById(): { [id: string]: string } | undefined {
  const { data, isLoading } = useCellTypes();

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
    url: `/api/ontology_tree_state?cellTypeId=%s`,
  },
  SOURCE_DATA: {
    queryKey: USE_SOURCE_DATA_QUERY,
    url: `/api/source_data?cellTypeId=%s`,
  },
  ENRICHED_GENES: {
    queryKey: USE_ENRICHED_GENES_QUERY,
    url: `/api/enriched_genes?cellTypeId=%s`,
  },
  CANONICAL_MARKERS: {
    queryKey: USE_CANONICAL_MARKERS_QUERY,
    url: `/api/canonical_markers?cellTypeId=%s`,
  },
  CL_DESCRIPTION: {
    queryKey: USE_CL_DESCRIPTION_QUERY,
    url: `/api/cl_description?cellTypeId=%s`,
  },
  DESCRIPTION: {
    queryKey: USE_DESCRIPTION_QUERY,
    url: `/api/description?cellTypeId=%s`,
  },
  CELL_CARDS: {
    queryKey: USE_CELL_CARDS_QUERY,
    url: "/api/cell_guides",
  },
};
