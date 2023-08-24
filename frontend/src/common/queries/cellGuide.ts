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
  CELL_GUIDE = "CELL_GUIDE",
  INITIAL_CELL_ONTOLOGY_TREE_TISSUE = "INITIAL_CELL_ONTOLOGY_TREE_TISSUE",
  TISSUE_CARDS = "TISSUE_CARDS",
  UBERON_DESCRIPTION = "UBERON_DESCRIPTION",
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
  | InitialCellOntologyTreeStateResponse
  | SourceDataQueryResponse
  | EnrichedGenesQueryResponse
  | CanonicalMarkersQueryResponse
  | ClDescriptionQueryResponse
  | DescriptionQueryResponse
  | CellGuideQueryResponse
  | UberonDescriptionQueryResponse;

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
export const USE_INITIAL_CELL_ONTOLOGY_TREE_STATE_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_INITIAL_CELL_ONTOLOGY_TREE_STATE],
  id: "cell-guide-cell-ontology-tree-state-query",
};

export interface TissueCountsPerCellType {
  [cell_type_id: string]: {
    n_cells: number;
    n_cells_rollup: number;
  };
}
export interface InitialCellOntologyTreeStateResponse {
  isExpandedNodes: string[];
  notShownWhenExpandedNodes: {
    [key: string]: string[];
  };
  tissueCounts?: TissueCountsPerCellType;
}

export const useCellOntologyTreeState = (
  entityId: string
): UseQueryResult<InitialCellOntologyTreeStateResponse> => {
  return useCellGuideQuery<InitialCellOntologyTreeStateResponse>(
    TYPES.INITIAL_CELL_ONTOLOGY_TREE,
    entityId
  );
};

/* ========== ontology_tree_state_tissue ========== */
export const USE_INITIAL_CELL_ONTOLOGY_TREE_STATE_TISSUE_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_INITIAL_CELL_ONTOLOGY_TREE_STATE_TISSUE],
  id: "cell-guide-cell-ontology-tree-state-tissue-query",
};

export const useCellOntologyTreeStateTissue = (
  entityId: string
): UseQueryResult<InitialCellOntologyTreeStateResponse> => {
  return useCellGuideQuery<InitialCellOntologyTreeStateResponse>(
    TYPES.INITIAL_CELL_ONTOLOGY_TREE_TISSUE,
    entityId
  );
};

/* ========== source_data ========== */
export const USE_SOURCE_DATA_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_SOURCE_DATA],
  id: "cell-guide-source-data-query",
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
  return useCellGuideQuery<SourceDataQueryResponse>(
    TYPES.SOURCE_DATA,
    entityId
  );
};

/* Tissues for cell type from source data dual mapping tissue name and ontology_term_id */
export function useSourceDataTissuesLabelToIdDualMap(
  cellTypeId: string
): Map<string, string> {
  const { data: sourceData } = useSourceData(cellTypeId);

  return useMemo(() => {
    if (!sourceData) return new Map<string, string>();

    const allTissuesMap = new Map<string, string>();
    sourceData.reduce((acc, curr) => {
      const { tissue } = curr;
      for (const t of tissue) {
        // map tissue label to ontology_term_id and ontology_term_id to tissue label
        acc.set(t.label, t.ontology_term_id);
        acc.set(t.ontology_term_id, t.label);
      }
      return acc;
    }, allTissuesMap);
    return allTissuesMap;
  }, [sourceData]);
}

/* ========== enriched_genes ========== */
export const USE_ENRICHED_GENES_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_ENRICHED_GENES],
  id: "cell-guide-enriched-genes-query",
};

interface EnrichedGenesQueryResponseEntry {
  me: number;
  pc: number;
  marker_score: number;
  symbol: string;
  name: string;
  organism: string;
  tissue: string;
}

export type EnrichedGenesQueryResponse = EnrichedGenesQueryResponseEntry[];

export const useEnrichedGenes = (
  entityId: string
): UseQueryResult<EnrichedGenesQueryResponse> => {
  return useCellGuideQuery<EnrichedGenesQueryResponse>(
    TYPES.ENRICHED_GENES,
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

/* ========== CL description ========== */
export const USE_CL_DESCRIPTION_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_CL_DESCRIPTION],
  id: "cell-guide-cl-description-query",
};

export type ClDescriptionQueryResponse = string;

export const useClDescription = (entityId: string): UseQueryResult<string> => {
  return useCellGuideQuery<string>(TYPES.CL_DESCRIPTION, entityId);
};

/* ========== UBERON description ========== */
export const USE_UBERON_DESCRIPTION_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_UBERON_DESCRIPTION],
  id: "cell-guide-uberon-description-query",
};

export type UberonDescriptionQueryResponse = string;

export const useUberonDescription = (
  entityId: string
): UseQueryResult<string> => {
  return useCellGuideQuery<string>(TYPES.UBERON_DESCRIPTION, entityId);
};

/* ========== description ========== */
export const USE_DESCRIPTION_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_DESCRIPTION],
  id: "cell-guide-description-query",
};

export type DescriptionQueryResponse = string;

export const useDescription = (entityId: string): UseQueryResult<string> => {
  return useCellGuideQuery<string>(TYPES.DESCRIPTION, entityId);
};

/* ========== cell_guide_cards ========== */
export const USE_CELL_GUIDE_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_CELL_GUIDE],
  id: "cell-guide-query",
};

interface CellGuideQueryResponseEntry {
  id: string;
  label: string;
  synonyms?: string[];
}

export type CellGuideQueryResponse = CellGuideQueryResponseEntry[];

export const useCellGuide = (): UseQueryResult<CellGuideQueryResponse> => {
  return useCellGuideQuery<CellGuideQueryResponse>(TYPES.CELL_GUIDE);
};

/* ========== cell types by Id ========== */
export function useCellTypesById():
  | { [cellTypeId: string]: CellGuideQueryResponseEntry }
  | undefined {
  const { data, isLoading } = useCellGuide();

  return useMemo(() => {
    if (!data || isLoading) return;

    const cellTypesById = {} as {
      [cellTypeId: string]: CellGuideQueryResponseEntry;
    };

    for (const cellType of data) {
      cellTypesById[cellType.id] = cellType;
    }

    return cellTypesById;
  }, [data, isLoading]);
}

/* ========== cell type names by Id ========== */

export function useCellTypeNamesById(): { [id: string]: string } | undefined {
  const { data, isLoading } = useCellGuide();

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
  entities: [ENTITIES.CELL_GUIDE_TISSUE_CARDS],
  id: "tissue-cards-query",
};

interface TissueCardsQueryResponseEntry {
  id: string;
  label: string;
}

export type TissueCardsQueryResponse = TissueCardsQueryResponseEntry[];

export const useTissueCards = (): UseQueryResult<TissueCardsQueryResponse> => {
  return useCellGuideQuery<TissueCardsQueryResponse>(TYPES.TISSUE_CARDS);
};

/* ========== all organs label to id map ========== */
export function useAllOrgansLabelToIdMap(): Map<string, string> {
  const { data: allOrgansData } = useTissueCards();
  return useMemo(() => {
    if (!allOrgansData) return new Map<string, string>();

    const allOrgansMap = new Map<string, string>();
    allOrgansData.reduce((acc, curr) => {
      const { id, label } = curr;
      acc.set(label, id);
      return acc;
    }, allOrgansMap);
    return allOrgansMap;
  }, [allOrgansData]);
}

/* ========== tissue descendants =========== */
/* export function useTissueDescendants(): Map<string, string[]> | undefined {
  return useMemo(() => {
    const tissueDescJsonObj = readJson(
      "src/components/common/Filter/descendant_mappings/tissue_descendants.json"
    );
    return new Map(Object.entries(tissueDescJsonObj));
  }, []);
}*/

/**
 * Mapping from data/response type to properties used for querying
 */
const QUERY_MAPPING: {
  [key in TYPES]: CellGuideQuery;
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
  UBERON_DESCRIPTION: {
    queryKey: USE_UBERON_DESCRIPTION_QUERY,
    url: `/api/uberon_description?entityId=%s`,
  },
  DESCRIPTION: {
    queryKey: USE_DESCRIPTION_QUERY,
    url: `/api/description?entityId=%s`,
  },
  CELL_GUIDE: {
    queryKey: USE_CELL_GUIDE_QUERY,
    url: "/api/cell_guide_cards",
  },
  TISSUE_CARDS: {
    queryKey: USE_TISSUE_CARDS_QUERY,
    url: "/api/tissue_cards",
  },
};
