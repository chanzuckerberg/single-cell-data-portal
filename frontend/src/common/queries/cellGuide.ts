import { useMemo, useEffect, useState } from "react";
import { useQuery, UseQueryResult } from "react-query";
import { DEFAULT_FETCH_OPTIONS, JSON_BODY_FETCH_OPTIONS } from "./common";
import { CELLGUIDE_DATA_URL } from "src/configs/configs";
import { ENTITIES } from "./entities";

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
  LATEST_SNAPSHOT_IDENTIFIER = "LATEST_SNAPSHOT_IDENTIFIER",
}

interface CellGuideQuery {
  queryKey: {
    entities: ENTITIES[];
    id: string;
  };
  urlSuffix: string;
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
  if (response.headers.get("Content-Length") === "0") {
    return undefined;
  }
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
  queryId = "", // Empty string if cell type is not needed for fetch function
  queryLatestSnapshotIdentifier = true
): UseQueryResult<T> {
  const { queryKey, urlSuffix } = QUERY_MAPPING[dataType];

  // if the query is "CL:0000000" make it "CL_0000000"
  queryId = queryId.replace(":", "_");

  const [latestSnapshotIdentifier, setLatestSnapshotIdentifier] = useState<
    LatestSnapshotIdentifierQueryResponse | undefined
  >(undefined);
  const {
    queryKey: queryKeyLatestSnapshotIdentifier,
    urlSuffix: urlSuffixLatestSnapshotIdentifier,
  } = QUERY_MAPPING[TYPES.LATEST_SNAPSHOT_IDENTIFIER];
  const { data: rawLatestSnapshotIdentifier } = useQuery(
    [queryKeyLatestSnapshotIdentifier],
    ({ signal }) =>
      fetchQuery({
        url: `${CELLGUIDE_DATA_URL}/${urlSuffixLatestSnapshotIdentifier.replace(
          "%s",
          queryId
        )}`,
        signal,
      }),
    {
      enabled: queryLatestSnapshotIdentifier,
      staleTime: Infinity,
    }
  );

  useEffect(() => {
    if (rawLatestSnapshotIdentifier) {
      setLatestSnapshotIdentifier(
        rawLatestSnapshotIdentifier as LatestSnapshotIdentifierQueryResponse
      );
    }
  }, [rawLatestSnapshotIdentifier]);

  const queryUrlSuffix = urlSuffix.replace("%s", queryId);
  const queryUrl = queryLatestSnapshotIdentifier
    ? `${CELLGUIDE_DATA_URL}/${latestSnapshotIdentifier}/${queryUrlSuffix}`
    : `${CELLGUIDE_DATA_URL}/${queryUrlSuffix}`;
  return useQuery(
    queryId ? [queryKey, queryId, latestSnapshotIdentifier] : [queryKey],
    ({ signal }) =>
      fetchQuery({
        url: queryUrl,
        signal,
      }),
    {
      enabled: !!latestSnapshotIdentifier,
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

export interface SourceCollectionsQueryResponseEntry {
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

export interface ComputationalMarkersQueryResponseEntry {
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

export const useComputationalMarkers = (
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

export interface CanonicalMarkersQueryResponseEntry {
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

/* ========== latest snapshot identifier ========== */
export const USE_LATEST_SNAPSHOT_IDENTIFIER_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_LATEST_SNAPSHOT_IDENTIFIER],
  id: "cell-guide-latest-snapshot-identifier-query",
};

export type LatestSnapshotIdentifierQueryResponse = string;

/* ========== description ========== */
export const USE_GPT_DESCRIPTION_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_GPT_DESCRIPTION],
  id: "cell-guide-gpt-description-query",
};

export type GptDescriptionQueryResponse = string;

export const useGptDescription = (entityId: string): UseQueryResult<string> => {
  return useCellGuideQuery<string>(TYPES.GPT_DESCRIPTION, entityId, false);
};

/* ========== SEO description ========== */
export const USE_GPT_SEO_DESCRIPTION_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_GPT_SEO_DESCRIPTION],
  id: "cell-guide-gpt-seo-description-query",
};

interface GptSeoDescriptionQueryResponse {
  name: string;
  description: string;
}

export const fetchGptSeoDescription = async (
  entityId: string
): Promise<GptSeoDescriptionQueryResponse> => {
  entityId = entityId.replace(":", "_");
  // This function is used server-side to fetch the GPT SEO description.
  const url = `${CELLGUIDE_DATA_URL}/${QUERY_MAPPING[
    TYPES.GPT_SEO_DESCRIPTION
  ].urlSuffix.replace("%s", entityId)}`;
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

interface CellTypeMetadataQueryResponse {
  [cellTypeId: string]: {
    id: string;
    name: string;
    synonyms?: string[];
    clDescription?: string;
  };
}

export const useCellTypeMetadata =
  (): UseQueryResult<CellTypeMetadataQueryResponse> => {
    return useCellGuideQuery<CellTypeMetadataQueryResponse>(
      TYPES.CELLTYPE_METADATA
    );
  };

/* ========== tissue_cards ========== */
export const USE_TISSUE_METADATA_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_TISSUE_METADATA],
  id: "cell-guide-tissue-metadata-query",
};

export interface TissueMetadataQueryResponse {
  [tissueId: string]: {
    id: string;
    name: string;
    synonyms?: string[];
    uberonDescription?: string;
  };
}

export const useTissueMetadata =
  (): UseQueryResult<TissueMetadataQueryResponse> => {
    return useCellGuideQuery<TissueMetadataQueryResponse>(
      TYPES.TISSUE_METADATA
    );
  };

export const fetchTissueMetadata =
  async (): Promise<TissueMetadataQueryResponse> => {
    // This function is used server-side to fetch the GPT SEO description.
    const latestSnapshotIdentifierUrl = `${CELLGUIDE_DATA_URL}/${
      QUERY_MAPPING[TYPES.LATEST_SNAPSHOT_IDENTIFIER].urlSuffix
    }`;
    const latestSnapshotIdentifierResponse = await fetch(
      latestSnapshotIdentifierUrl
    );
    const latestSnapshotIdentifier =
      await latestSnapshotIdentifierResponse.text();
    const url = `${CELLGUIDE_DATA_URL}/${latestSnapshotIdentifier}/${
      QUERY_MAPPING[TYPES.TISSUE_METADATA].urlSuffix
    }`;
    const response = await fetch(url);
    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }
    return await response.json();
  };

/* ========== Lookup tables for organs ========== */
export function useAllOrgansLookupTables(): {
  organLabelToIdMap: Map<string, string>;
} {
  const { data: allOrgansData } = useTissueMetadata();
  return useMemo(() => {
    if (!allOrgansData) {
      return {
        organLabelToIdMap: new Map<string, string>(),
      };
    }

    const allOrgansLabelToIdMap = new Map<string, string>();
    for (const organId in allOrgansData) {
      const organData = allOrgansData[organId];
      allOrgansLabelToIdMap.set(organData.name, organData.id);
    }
    return {
      organLabelToIdMap: allOrgansLabelToIdMap,
    };
  }, [allOrgansData]);
}

/* ========== Lookup tables for tissues ========== */
export function useAllTissuesLookupTables(cellTypeId: string): {
  allTissuesLabelToIdMap: Map<string, string>;
} {
  const { data: sourceData } = useSourceData(cellTypeId);

  return useMemo(() => {
    if (!sourceData) {
      return { allTissuesLabelToIdMap: new Map<string, string>() };
    }

    const allTissuesLabelToIdLookup = new Map<string, string>();

    for (const source of sourceData) {
      const tissueList = source.tissue;
      for (const tissue of tissueList) {
        allTissuesLabelToIdLookup.set(tissue.label, tissue.ontology_term_id);
      }
    }
    return { allTissuesLabelToIdMap: allTissuesLabelToIdLookup };
  }, [sourceData]);
}

/**
 * Mapping from data/response type to properties used for querying
 */
const QUERY_MAPPING: {
  [key in TYPES]: CellGuideQuery;
} = {
  CELL_ONTOLOGY_TREE: {
    queryKey: USE_CELL_ONTOLOGY_TREE_QUERY,
    urlSuffix: "ontology_graph.json",
  },
  CELL_ONTOLOGY_TREE_STATE_CELLTYPE: {
    queryKey: USE_CELL_ONTOLOGY_TREE_STATE_CELLTYPE_QUERY,
    urlSuffix: `cell_type_ontology_tree_state/%s.json`,
  },
  CELL_ONTOLOGY_TREE_STATE_TISSUE: {
    queryKey: USE_CELL_ONTOLOGY_TREE_STATE_TISSUE_QUERY,
    urlSuffix: `tissue_ontology_tree_state/%s.json`,
  },
  SOURCE_COLLECTIONS: {
    queryKey: USE_SOURCE_COLLECTIONS_QUERY,
    urlSuffix: `source_collections/%s.json`,
  },
  COMPUTATIONAL_MARKERS: {
    queryKey: USE_COMPUTATIONAL_MARKERS_QUERY,
    urlSuffix: `computational_marker_genes/%s.json`,
  },
  CANONICAL_MARKERS: {
    queryKey: USE_CANONICAL_MARKERS_QUERY,
    urlSuffix: `canonical_marker_genes/%s.json`,
  },
  GPT_DESCRIPTION: {
    queryKey: USE_GPT_DESCRIPTION_QUERY,
    urlSuffix: `gpt_descriptions/%s.json`,
  },
  GPT_SEO_DESCRIPTION: {
    queryKey: USE_GPT_SEO_DESCRIPTION_QUERY,
    urlSuffix: `gpt_seo_descriptions/%s.json`,
  },
  CELLTYPE_METADATA: {
    queryKey: USE_CELLTYPE_METADATA_QUERY,
    urlSuffix: "celltype_metadata.json",
  },
  TISSUE_METADATA: {
    queryKey: USE_TISSUE_METADATA_QUERY,
    urlSuffix: "tissue_metadata.json",
  },
  LATEST_SNAPSHOT_IDENTIFIER: {
    queryKey: USE_LATEST_SNAPSHOT_IDENTIFIER_QUERY,
    urlSuffix: "latest_snapshot_identifier",
  },
};
