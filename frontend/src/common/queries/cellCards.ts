import { useMemo } from "react";
import { useQuery, UseQueryResult } from "react-query";
import { DEFAULT_FETCH_OPTIONS, JSON_BODY_FETCH_OPTIONS } from "./common";
import { ENTITIES } from "./entities";

// source_data
interface WikipediaDescriptionQueryResponse {
  content: string;
}

async function fetchWikipediaDescriptionQuery({
  cellTypeId,
  signal,
}: {
  cellTypeId: string;
  signal?: AbortSignal;
}): Promise<WikipediaDescriptionQueryResponse | undefined> {
  const url = `/api/scrape?cellTypeId=${cellTypeId}`;
  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    method: "GET",
    signal,
  });
  if (response.status === 404) return undefined;
  const json: WikipediaDescriptionQueryResponse = await response.json();

  if (!response.ok) {
    throw json;
  }

  return json;
}

export const USE_WIKIPEDIA_DESCRIPTION_QUERY = {
  entities: [ENTITIES.CELL_EXPLORER_WIKIPEDIA_DESCRIPTION],
  id: "cell-explorer-wikipedia-description-query",
};

export function useWikipediaDescription(
  cellTypeId: string
): UseQueryResult<WikipediaDescriptionQueryResponse> {
  return useQuery(
    [USE_WIKIPEDIA_DESCRIPTION_QUERY, cellTypeId],
    ({ signal }) => fetchWikipediaDescriptionQuery({ cellTypeId, signal }),
    {
      enabled: true,
      staleTime: Infinity,
    }
  );
}

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
  if (response.status === 404) return undefined;
  const json: SourceDataQueryResponse = await response.json();

  if (!response.ok) {
    throw json;
  }

  return json;
}

export const USE_SOURCE_DATA_QUERY = {
  entities: [ENTITIES.CELL_EXPLORER_SOURCE_DATA],
  id: "cell-explorer-source-data-query",
};

export function useSourceData(
  cellTypeId: string
): UseQueryResult<SourceDataQueryResponse> {
  return useQuery(
    [USE_SOURCE_DATA_QUERY, cellTypeId],
    ({ signal }) => fetchSourceDataQuery({ cellTypeId, signal }),
    {
      enabled: true,
      staleTime: Infinity,
    }
  );
}

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
  if (response.status === 404) return undefined;
  const json: EnrichedGenesQueryResponse = await response.json();

  if (!response.ok) {
    throw json;
  }

  return json;
}

export const USE_ENRICHED_GENES_QUERY = {
  entities: [ENTITIES.CELL_EXPLORER_ENRICHED_GENES],
  id: "cell-explorer-enriched-genes-query",
};

export function useEnrichedGenes(
  cellTypeId: string
): UseQueryResult<EnrichedGenesQueryResponse> {
  return useQuery(
    [USE_ENRICHED_GENES_QUERY, cellTypeId],
    ({ signal }) => fetchEnrichedGenesQuery({ cellTypeId, signal }),
    {
      enabled: true,
      staleTime: Infinity,
    }
  );
}

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
  if (response.status === 404) return undefined;
  const json: CanonicalMarkersQueryResponse = await response.json();

  if (!response.ok) {
    throw json;
  }

  return json;
}

export const USE_CANONICAL_MARKERS_QUERY = {
  entities: [ENTITIES.CELL_EXPLORER_CANONICAL_MARKERS],
  id: "cell-explorer-canonical-markersquery",
};

export function useCanonicalMarkers(
  cellTypeId: string
): UseQueryResult<CanonicalMarkersQueryResponse> {
  return useQuery(
    [USE_CANONICAL_MARKERS_QUERY, cellTypeId],
    ({ signal }) => fetchCanonicalMarkersQuery({ cellTypeId, signal }),
    {
      enabled: true,
      staleTime: Infinity,
    }
  );
}

// description
async function fetchDescription({
  cellTypeId,
  signal,
}: {
  cellTypeId: string;
  signal?: AbortSignal;
}): Promise<string | undefined> {
  const url = `/api/description?cellTypeId=${cellTypeId}`;
  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    method: "GET",
    signal,
  });
  if (response.status === 404) return undefined;
  const json: string = await response.json();

  if (!response.ok) {
    throw json;
  }
  return json;
}

export const USE_DESCRIPTION_QUERY = {
  entities: [ENTITIES.CELL_EXPLORER_DESCRIPTION],
  id: "cell-explorer-description-query",
};

export function useDescription(cellTypeId: string): UseQueryResult<string> {
  return useQuery(
    [USE_DESCRIPTION_QUERY, cellTypeId],
    ({ signal }) => fetchDescription({ cellTypeId, signal }),
    {
      enabled: true,
      staleTime: Infinity,
    }
  );
}

// cell_cards
interface CellCardsQueryResponseEntry {
  id: string;
  label: string;
}

type CellCardsQueryResponse = CellCardsQueryResponseEntry[];

async function fetchCellCardsQuery(
  signal?: AbortSignal
): Promise<CellCardsQueryResponse | undefined> {
  const url = "/api/cell_cards";
  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    method: "GET",
    signal,
  });
  if (response.status === 404) return undefined;
  const json: CellCardsQueryResponse = await response.json();

  if (!response.ok) {
    throw json;
  }

  return json;
}

export const USE_CELL_CARDS_QUERY = {
  entities: [ENTITIES.CELL_EXPLORER_CELL_CARDS],
  id: "cell-cards-query",
};

export function useCellTypes(): UseQueryResult<CellCardsQueryResponse> {
  return useQuery(
    [USE_CELL_CARDS_QUERY],
    ({ signal }) => fetchCellCardsQuery(signal),
    {
      enabled: true,
      staleTime: Infinity,
    }
  );
}

export function useCellTypesById() {
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
