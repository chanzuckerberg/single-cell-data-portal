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
  queryFn: ({
    cellTypeId,
    signal,
  }: {
    cellTypeId: string;
    signal: AbortSignal | undefined;
  }) => Promise<CELL_CARD_RESPONSE_TYPE | undefined>;
}

type CELL_CARD_RESPONSE_TYPE =
  | CellOntologyTreeResponse
  | InitialCellOntologyTreeStateResponse
  | SourceDataQueryResponse
  | EnrichedGenesQueryResponse
  | CanonicalMarkersQueryResponse
  | ClDescriptionQueryResponse
  | DescriptionQueryResponse
  | CellCardsQueryResponse;

export function useCellCardQuery(
  dataType: TYPES,
  cellTypeId = "" // Empty string if cell type is not needed for fetch function
): UseQueryResult<CELL_CARD_RESPONSE_TYPE> {
  const { queryKey, queryFn } = QUERY_MAPPING[dataType];

  return useQuery(
    cellTypeId ? [queryKey, cellTypeId] : [queryKey],
    ({ signal }) => queryFn({ cellTypeId, signal }),
    {
      enabled: true,
      staleTime: Infinity,
    }
  );
}

// ontology_tree
export interface CellOntologyTreeResponse {
  name: string;
  id: string;
  n_cells_rollup: number;
  n_cells: number;
  children?: this[];
  _children?: this[];
}

async function fetchOntologyTreeQuery({
  signal,
}: {
  cellTypeId?: string;
  signal?: AbortSignal;
}): Promise<CellOntologyTreeResponse | undefined> {
  const url = `/api/ontology_tree`;
  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    method: "GET",
    signal,
  });
  if (response.status === 204) return undefined;
  const json: CellOntologyTreeResponse = await response.json();

  if (!response.ok) {
    throw json;
  }

  return json;
}

export const USE_CELL_ONTOLOGY_TREE_QUERY = {
  entities: [ENTITIES.CELL_CARDS_CELL_ONTOLOGY_TREE],
  id: "cell-explorer-cell-ontology-tree-query",
};

// ontology_tree_state
export interface InitialCellOntologyTreeStateResponse {
  isExpandedNodes: string[];
  notShownWhenExpandedNodes: {
    [key: string]: string[];
  };
}

async function fetchOntologyTreeStateQuery({
  cellTypeId,
  signal,
}: {
  cellTypeId: string;
  signal?: AbortSignal;
}): Promise<InitialCellOntologyTreeStateResponse | undefined> {
  const url = `/api/ontology_tree_state?cellTypeId=${cellTypeId}`;
  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    method: "GET",
    signal,
  });
  if (response.status === 204) return undefined;
  const json: InitialCellOntologyTreeStateResponse = await response.json();

  if (!response.ok) {
    throw json;
  }

  return json;
}

export const USE_INITIAL_CELL_ONTOLOGY_TREE_STATE_QUERY = {
  entities: [ENTITIES.CELL_CARDS_INITIAL_CELL_ONTOLOGY_TREE_STATE],
  id: "cell-explorer-cell-ontology-tree-query",
};

// source_data
interface SourceDataQueryResponseEntry {
  collection_name: string;
  collection_url: string;
  publication_url: string;
  publication_title: string;
  tissue: { label: string; ontology_term_id: string }[];
  disease: { label: string; ontology_term_id: string }[];
  organism: { label: string; ontology_term_id: string }[];
}

type SourceDataQueryResponse = SourceDataQueryResponseEntry[];

async function fetchSourceDataQuery({
  cellTypeId,
  signal,
}: {
  cellTypeId: string;
  signal?: AbortSignal;
}): Promise<SourceDataQueryResponse | undefined> {
  const url = `/api/source_data?cellTypeId=${cellTypeId}`;
  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    method: "GET",
    signal,
  });
  if (response.status === 204) return undefined;
  const json: SourceDataQueryResponse = await response.json();

  if (!response.ok) {
    throw json;
  }

  return json;
}

export const USE_SOURCE_DATA_QUERY = {
  entities: [ENTITIES.CELL_CARDS_SOURCE_DATA],
  id: "cell-explorer-source-data-query",
};

// enriched_genes
interface EnrichedGenesQueryResponseEntry {
  me: number;
  pc: number;
  symbol: string;
  name: string;
  organism: string;
}

type EnrichedGenesQueryResponse = EnrichedGenesQueryResponseEntry[];

async function fetchEnrichedGenesQuery({
  cellTypeId,
  signal,
}: {
  cellTypeId: string;
  signal?: AbortSignal;
}): Promise<EnrichedGenesQueryResponse | undefined> {
  const url = `/api/enriched_genes?cellTypeId=${cellTypeId}`;
  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    method: "GET",
    signal,
  });
  if (response.status === 204) return undefined;
  const json: EnrichedGenesQueryResponse = await response.json();

  if (!response.ok) {
    throw json;
  }

  return json;
}

export const USE_ENRICHED_GENES_QUERY = {
  entities: [ENTITIES.CELL_CARDS_ENRICHED_GENES],
  id: "cell-explorer-enriched-genes-query",
};

// canonical_markers
interface CanonicalMarkersQueryResponseEntry {
  tissue_general: string;
  tissue_specific: string;
  symbol: string;
  name: string;
  publication: string;
  publication_titles: string;
}

type CanonicalMarkersQueryResponse = CanonicalMarkersQueryResponseEntry[];

async function fetchCanonicalMarkersQuery({
  cellTypeId,
  signal,
}: {
  cellTypeId: string;
  signal?: AbortSignal;
}): Promise<CanonicalMarkersQueryResponse | undefined> {
  const url = `/api/canonical_markers?cellTypeId=${cellTypeId}`;
  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    method: "GET",
    signal,
  });
  if (response.status === 204) return undefined;
  const json: CanonicalMarkersQueryResponse = await response.json();

  if (!response.ok) {
    throw json;
  }

  return json;
}

export const USE_CANONICAL_MARKERS_QUERY = {
  entities: [ENTITIES.CELL_CARDS_CANONICAL_MARKERS],
  id: "cell-explorer-canonical-markersquery",
};

// CL description
async function fetchClDescription({
  cellTypeId,
  signal,
}: {
  cellTypeId: string;
  signal?: AbortSignal;
}): Promise<ClDescriptionQueryResponse | undefined> {
  const url = `/api/cl_description?cellTypeId=${cellTypeId}`;
  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    method: "GET",
    signal,
  });
  if (response.status === 204) return undefined;
  const json: string = await response.json();

  if (!response.ok) {
    throw json;
  }
  return json;
}

export const USE_CL_DESCRIPTION_QUERY = {
  entities: [ENTITIES.CELL_CARDS_CL_DESCRIPTION],
  id: "cell-explorer-cl-description-query",
};

type ClDescriptionQueryResponse = string;

// description
async function fetchDescription({
  cellTypeId,
  signal,
}: {
  cellTypeId: string;
  signal?: AbortSignal;
}): Promise<DescriptionQueryResponse | undefined> {
  const url = `/api/description?cellTypeId=${cellTypeId}`;
  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    method: "GET",
    signal,
  });
  if (response.status === 204) return undefined;
  const json: string = await response.json();

  if (!response.ok) {
    throw json;
  }
  return json;
}

export const USE_DESCRIPTION_QUERY = {
  entities: [ENTITIES.CELL_CARDS_DESCRIPTION],
  id: "cell-explorer-description-query",
};

type DescriptionQueryResponse = string;

// cell_guides
interface CellCardsQueryResponseEntry {
  id: string;
  label: string;
}

type CellCardsQueryResponse = CellCardsQueryResponseEntry[];

async function fetchCellCardsQuery({
  signal,
}: {
  cellTypeId?: string;
  signal?: AbortSignal;
}): Promise<CellCardsQueryResponse | undefined> {
  const url = "/api/cell_guides";
  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    method: "GET",
    signal,
  });
  if (response.status === 204) return undefined;
  const json: CellCardsQueryResponse = await response.json();

  if (!response.ok) {
    throw json;
  }

  return json;
}

export const USE_CELL_CARDS_QUERY = {
  entities: [ENTITIES.CELL_CARDS_CELL_CARDS],
  id: "cell-cards-query",
};

export function useCellTypesById(): { [id: string]: string } | undefined {
  const { data, isLoading } = useCellCardQuery(
    TYPES.CELL_CARDS
  ) as UseQueryResult<CellCardsQueryResponse>; //useCellTypes();

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

const QUERY_MAPPING: {
  [key in TYPES]: CellCardQuery;
} = {
  CELL_ONTOLOGY_TREE: {
    queryKey: USE_CELL_ONTOLOGY_TREE_QUERY,
    queryFn: fetchOntologyTreeQuery,
  },
  INITIAL_CELL_ONTOLOGY_TREE: {
    queryKey: USE_INITIAL_CELL_ONTOLOGY_TREE_STATE_QUERY,
    queryFn: fetchOntologyTreeStateQuery,
  },
  SOURCE_DATA: {
    queryKey: USE_SOURCE_DATA_QUERY,
    queryFn: fetchSourceDataQuery,
  },
  ENRICHED_GENES: {
    queryKey: USE_ENRICHED_GENES_QUERY,
    queryFn: fetchEnrichedGenesQuery,
  },
  CANONICAL_MARKERS: {
    queryKey: USE_CANONICAL_MARKERS_QUERY,
    queryFn: fetchCanonicalMarkersQuery,
  },
  CL_DESCRIPTION: {
    queryKey: USE_CL_DESCRIPTION_QUERY,
    queryFn: fetchClDescription,
  },
  DESCRIPTION: {
    queryKey: USE_DESCRIPTION_QUERY,
    queryFn: fetchDescription,
  },
  CELL_CARDS: {
    queryKey: USE_CELL_CARDS_QUERY,
    queryFn: fetchCellCardsQuery,
  },
};
