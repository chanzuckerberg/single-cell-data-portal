import { useMemo } from "react";
import { useQuery, UseQueryResult } from "react-query";
import { DEFAULT_FETCH_OPTIONS, JSON_BODY_FETCH_OPTIONS } from "./common";
import { ENTITIES } from "./entities";

const CELLGUIDE_DATA_URL_PREFIX = "cellguide";

function getCellGuideDataUrl(urlSuffix: string) {
  const currentUrl = new URL(window.location.href);
  const [, ...domainParts] = currentUrl.hostname.split(".");
  const baseDomain = domainParts.join(".");
  currentUrl.hostname = `${CELLGUIDE_DATA_URL_PREFIX}.${baseDomain}`;
  return `${currentUrl.toString()}${urlSuffix}`;
}

export enum TYPES {
  CELL_ONTOLOGY_TREE = "CELL_ONTOLOGY_TREE",
  CELL_ONTOLOGY_TREE_STATE_CELLTYPE = "CELL_ONTOLOGY_TREE_STATE_CELLTYPE",
  SOURCE_COLLECTIONS = "SOURCE_COLLECTIONS",
  COMPUTATIONAL_MARKERS = "COMPUTATIONAL_MARKERS",
  CANONICAL_MARKERS = "CANONICAL_MARKERS",
  GPT_DESCRIPTION = "GPT_DESCRIPTION",
  CELL_ONTOLOGY_TREE_STATE_TISSUE = "CELL_ONTOLOGY_TREE_STATE_TISSUE",
  TISSUE_METADATA = "TISSUE_METADATA",
  CELLTYPE_METADATA = "CELLTYPE_METADATA",
  GPT_SEO_DESCRIPTION = "GPT_SEO_DESCRIPTION",
}

interface CellGuideQuery {
  queryKey: {
    entities: ENTITIES[];
    id: string;
  };
  url: string;
}

export type CellGuideResponse =
  | CellOntologyTreeResponse
  | CellOntologyTreeStateResponse
  | SourceCollectionsQueryResponse
  | ComputationalMarkersQueryResponse
  | CanonicalMarkersQueryResponse
  | GptDescriptionQueryResponse
  | CellTypeMetadataQueryResponse
  | TissueMetadataQueryResponse;

/**
 * Generic fetch function
 */
async function fetchQuery({
  url,
  signal,
}: {
  url: string;
  signal: AbortSignal | undefined;
}): Promise<CellGuideResponse | undefined> {
  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    method: "GET",
    signal,
  });
  if (response.status === 204) return undefined;
  const json: CellGuideResponse = await response.json();

  if (!response.ok) {
    throw json;
  }

  return json;
}

/**
 * Generic cell guide hook
 */
export function useCellGuideQuery<T = CellGuideResponse>(
  dataType: TYPES,
  queryId = "" // Empty string if cell type is not needed for fetch function
): UseQueryResult<T> {
  const { queryKey, url: rawUrl } = QUERY_MAPPING[dataType];
  return useQuery(
    queryId ? [queryKey, queryId] : [queryKey],
    ({ signal }) =>
      fetchQuery({
        url: getCellGuideDataUrl(rawUrl.replace("%s", queryId)), // Replacing raw url with entityId if applicable
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
  entities: [ENTITIES.CELL_GUIDE_CELL_ONTOLOGY_TREE],
  id: "cell-guide-cell-ontology-tree-query",
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
    return useCellGuideQuery<CellOntologyTreeResponse>(
      TYPES.CELL_ONTOLOGY_TREE
    );
  };

/* ========== ontology_tree_state ========== */
export const USE_CELL_ONTOLOGY_TREE_STATE_CELLTYPE_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_CELL_ONTOLOGY_TREE_STATE_CELLTYPE],
  id: "cell-guide-cell-ontology-tree-state-celltype-query",
};

export interface TissueCountsPerCellType {
  [cell_type_id: string]: {
    n_cells: number;
    n_cells_rollup: number;
  };
}
export interface CellOntologyTreeStateResponse {
  isExpandedNodes: string[];
  notShownWhenExpandedNodes: {
    [key: string]: string[];
  };
  tissueCounts?: TissueCountsPerCellType;
}

export const useCellOntologyTreeStateCellType = (
  entityId: string
): UseQueryResult<CellOntologyTreeStateResponse> => {
  return useCellGuideQuery<CellOntologyTreeStateResponse>(
    TYPES.CELL_ONTOLOGY_TREE_STATE_CELLTYPE,
    entityId
  );
};

/* ========== ontology_tree_state_tissue ========== */
export const USE_CELL_ONTOLOGY_TREE_STATE_TISSUE_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_CELL_ONTOLOGY_TREE_STATE_TISSUE],
  id: "cell-guide-cell-ontology-tree-state-tissue-query",
};

export const useCellOntologyTreeStateTissue = (
  entityId: string
): UseQueryResult<CellOntologyTreeStateResponse> => {
  return useCellGuideQuery<CellOntologyTreeStateResponse>(
    TYPES.CELL_ONTOLOGY_TREE_STATE_TISSUE,
    entityId
  );
};

/* ========== source_data ========== */
export const USE_SOURCE_COLLECTIONS_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_SOURCE_COLLECTIONS],
  id: "cell-guide-source-collections-query",
};

interface SourceCollectionsQueryResponseEntry {
  collection_name: string;
  collection_url: string;
  publication_url: string;
  publication_title: string;
  tissue: { label: string; ontology_term_id: string }[];
  disease: { label: string; ontology_term_id: string }[];
  organism: { label: string; ontology_term_id: string }[];
}

export type SourceCollectionsQueryResponse =
  SourceCollectionsQueryResponseEntry[];

export const useSourceData = (
  entityId: string
): UseQueryResult<SourceCollectionsQueryResponse> => {
  return useCellGuideQuery<SourceCollectionsQueryResponse>(
    TYPES.SOURCE_COLLECTIONS,
    entityId
  );
};

/* ========== enriched_genes ========== */
export const USE_COMPUTATIONAL_MARKERS_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_COMPUTATIONAL_MARKERS],
  id: "cell-guide-computational-markers-query",
};

interface ComputationalMarkersQueryResponseEntry {
  me: number;
  pc: number;
  marker_score: number;
  symbol: string;
  name: string;
  groupby_dims: {
    organism_ontology_term_label: string;
    tissue_ontology_term_label?: string;
  };
}

export type ComputationalMarkersQueryResponse =
  ComputationalMarkersQueryResponseEntry[];

export const useEnrichedGenes = (
  entityId: string
): UseQueryResult<ComputationalMarkersQueryResponse> => {
  return useCellGuideQuery<ComputationalMarkersQueryResponse>(
    TYPES.COMPUTATIONAL_MARKERS,
    entityId
  );
};

/* ========== canonical_markers ========== */
export const USE_CANONICAL_MARKERS_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_CANONICAL_MARKERS],
  id: "cell-guide-canonical-markersquery",
};

interface CanonicalMarkersQueryResponseEntry {
  tissue: string;
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
  return useCellGuideQuery<CanonicalMarkersQueryResponse>(
    TYPES.CANONICAL_MARKERS,
    entityId
  );
};

/* ========== description ========== */
export const USE_GPT_DESCRIPTION_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_GPT_DESCRIPTION],
  id: "cell-guide-gpt-description-query",
};

export type GptDescriptionQueryResponse = string;

export const useGptDescription = (entityId: string): UseQueryResult<string> => {
  return useCellGuideQuery<string>(TYPES.GPT_DESCRIPTION, entityId);
};

/* ========== SEO description ========== */
export const USE_GPT_SEO_DESCRIPTION_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_GPT_SEO_DESCRIPTION],
  id: "cell-guide-gpt-seo-description-query",
};

interface GptSeoDescriptionQueryResponse {
  [cellTypeId: string]: { name: string; description: string };
}

export const fetchGptSeoDescription = async (
  entityId: string
): Promise<GptSeoDescriptionQueryResponse> => {
  // This function is used server-side to fetch the GPT SEO description.
  const url = getCellGuideDataUrl(
    QUERY_MAPPING[TYPES.GPT_SEO_DESCRIPTION].url.replace("%s", entityId)
  );
  const response = await fetch(url);
  if (!response.ok) {
    throw new Error(`HTTP error! status: ${response.status}`);
  }
  return await response.json();
};

/* ========== cell_guide_cards ========== */
export const USE_CELLTYPE_METADATA_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_CELLTYPE_METADATA],
  id: "cell-guide-celltype-metadata-query",
};

interface CellTypeMetadataQueryResponseEntry {
  id: string;
  label: string;
  synonyms?: string[];
  clDescription?: string;
}

export type CellTypeMetadataQueryResponse =
  CellTypeMetadataQueryResponseEntry[];

export const useCellTypeMetadata =
  (): UseQueryResult<CellTypeMetadataQueryResponse> => {
    return useCellGuideQuery<CellTypeMetadataQueryResponse>(
      TYPES.CELLTYPE_METADATA
    );
  };

/* ========== cell types by Id ========== */
export function useCellTypesById():
  | { [cellTypeId: string]: CellTypeMetadataQueryResponseEntry }
  | undefined {
  const { data, isLoading } = useCellTypeMetadata();

  return useMemo(() => {
    if (!data || isLoading) return;

    const cellTypesById = {} as {
      [cellTypeId: string]: CellTypeMetadataQueryResponseEntry;
    };

    for (const cellType of data) {
      cellTypesById[cellType.id] = cellType;
    }

    return cellTypesById;
  }, [data, isLoading]);
}

/* ========== cell type names by Id ========== */

export function useCellTypeNamesById(): { [id: string]: string } | undefined {
  const { data, isLoading } = useCellTypeMetadata();

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
export const USE_TISSUE_METADATA_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_TISSUE_METADATA],
  id: "cell-guide-tissue-metadata-query",
};

interface TissueMetadataQueryResponseEntry {
  id: string;
  label: string;
  synonyms?: string[];
  uberonDescription?: string;
}

export type TissueMetadataQueryResponse = TissueMetadataQueryResponseEntry[];

export const useTissueMetadata =
  (): UseQueryResult<TissueMetadataQueryResponse> => {
    return useCellGuideQuery<TissueMetadataQueryResponse>(
      TYPES.TISSUE_METADATA
    );
  };

export const fetchTissueMetadata =
  async (): Promise<TissueMetadataQueryResponse> => {
    // This function is used server-side to fetch the GPT SEO description.
    const url = getCellGuideDataUrl(QUERY_MAPPING[TYPES.TISSUE_METADATA].url);
    const response = await fetch(url);
    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }
    return await response.json();
  };

/**
 * Mapping from data/response type to properties used for querying
 */
const QUERY_MAPPING: {
  [key in TYPES]: CellGuideQuery;
} = {
  CELL_ONTOLOGY_TREE: {
    queryKey: USE_CELL_ONTOLOGY_TREE_QUERY,
    url: "/ontology_graph.json",
  },
  CELL_ONTOLOGY_TREE_STATE_CELLTYPE: {
    queryKey: USE_CELL_ONTOLOGY_TREE_STATE_CELLTYPE_QUERY,
    url: `/cell_type_ontology_tree_state/%s`,
  },
  CELL_ONTOLOGY_TREE_STATE_TISSUE: {
    queryKey: USE_CELL_ONTOLOGY_TREE_STATE_TISSUE_QUERY,
    url: `/tissue_ontology_tree_state/%s`,
  },
  SOURCE_COLLECTIONS: {
    queryKey: USE_SOURCE_COLLECTIONS_QUERY,
    url: `/source_collections/%s`,
  },
  COMPUTATIONAL_MARKERS: {
    queryKey: USE_COMPUTATIONAL_MARKERS_QUERY,
    url: `/computational_marker_genes/%s`,
  },
  CANONICAL_MARKERS: {
    queryKey: USE_CANONICAL_MARKERS_QUERY,
    url: `/canonical_marker_genes/%s`,
  },
  GPT_DESCRIPTION: {
    queryKey: USE_GPT_DESCRIPTION_QUERY,
    url: `/gpt_description/%s`,
  },
  GPT_SEO_DESCRIPTION: {
    queryKey: USE_GPT_SEO_DESCRIPTION_QUERY,
    url: `/gpt_description/%s`,
  },
  CELLTYPE_METADATA: {
    queryKey: USE_CELLTYPE_METADATA_QUERY,
    url: "/celltype_metadata.json",
  },
  TISSUE_METADATA: {
    queryKey: USE_TISSUE_METADATA_QUERY,
    url: "/tissue_metadata.json",
  },
};
