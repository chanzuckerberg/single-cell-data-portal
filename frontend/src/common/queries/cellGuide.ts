import pako from "pako";
import { useMemo, useEffect, useState } from "react";
import { useQuery, UseQueryResult } from "react-query";
import { DEFAULT_FETCH_OPTIONS, JSON_BODY_FETCH_OPTIONS } from "./common";
import { CELLGUIDE_DATA_URL, API_URL } from "src/configs/configs";
import { ENTITIES } from "./entities";
import { EMPTY_ARRAY } from "../constants/utils";

export enum TYPES {
  CELL_ONTOLOGY_TREE = "CELL_ONTOLOGY_TREE",
  CELL_ONTOLOGY_TREE_STATE_CELLTYPE = "CELL_ONTOLOGY_TREE_STATE_CELLTYPE",
  SOURCE_COLLECTIONS = "SOURCE_COLLECTIONS",
  COMPUTATIONAL_MARKERS = "COMPUTATIONAL_MARKERS",
  CANONICAL_MARKERS = "CANONICAL_MARKERS",
  GPT_DESCRIPTION = "GPT_DESCRIPTION",
  VALIDATED_DESCRIPTION = "VALIDATED_DESCRIPTION",
  CELL_ONTOLOGY_TREE_STATE_TISSUE = "CELL_ONTOLOGY_TREE_STATE_TISSUE",
  TISSUE_METADATA = "TISSUE_METADATA",
  CELLTYPE_METADATA = "CELLTYPE_METADATA",
  MARKER_GENE_PRESENCE = "MARKER_GENE_PRESENCE",
  GPT_SEO_DESCRIPTION = "GPT_SEO_DESCRIPTION",
  CELLTYPE_TISSUE_MAPPING = "CELLTYPE_TISSUE_MAPPING",
  LATEST_SNAPSHOT_IDENTIFIER = "LATEST_SNAPSHOT_IDENTIFIER",
  VALID_EXPLORER_CXGS = "VALID_EXPLORER_CXGS",
}

// Suffix the cellguide data url with the remote dev prefix
const IS_RDEV = API_URL.includes(".rdev.single-cell.czi.technology");
let CELLGUIDE_DATA_URL_WITH_RDEV_SUFFIX = CELLGUIDE_DATA_URL;
if (IS_RDEV) {
  const REMOTE_DEV_PREFIX = API_URL.split("//")[1].split("-backend")[0];
  CELLGUIDE_DATA_URL_WITH_RDEV_SUFFIX = `${CELLGUIDE_DATA_URL}/env-rdev-cellguide/${REMOTE_DEV_PREFIX}`;
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

  // need to include access-controle-expose-headers in the response header
  // to use Content-Encoding. For now, check if file ends with .gz.
  if (url.endsWith(".gz")) {
    const arrayBuffer = await response.arrayBuffer();
    const decompressedData = pako.inflate(new Uint8Array(arrayBuffer), {
      to: "string",
    });
    const json: CellGuideResponse = JSON.parse(decompressedData);
    return json;
  } else {
    const json: CellGuideResponse = await response.json();
    return json;
  }
}

/**
 * Generic cell guide hook
 */
// Mapping from organism name to taxon id. Hardcode this for now.
// This will need to be updated if Census ever includes more organisms.
export const ORGANISM_NAME_TO_TAXON_ID_MAPPING = {
  "Homo sapiens": "NCBITaxon_9606",
  "Mus musculus": "NCBITaxon_10090",
};

interface CellGuideQueryProps {
  dataType: TYPES;
  queryId?: string;
  organismName?: string;
  queryLatestSnapshotIdentifier?: boolean;
}
export function useCellGuideQuery<T = CellGuideResponse>({
  dataType,
  queryId = "",
  organismName = "",
  queryLatestSnapshotIdentifier = true,
}: CellGuideQueryProps): UseQueryResult<T> {
  // if organismName is a key in ORGANISM_NAME_TO_TAXON_ID_MAPPING, set it to the value
  const organismId =
    ORGANISM_NAME_TO_TAXON_ID_MAPPING[
      organismName as keyof typeof ORGANISM_NAME_TO_TAXON_ID_MAPPING
    ] || organismName;

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
        url: `${CELLGUIDE_DATA_URL_WITH_RDEV_SUFFIX}/${urlSuffixLatestSnapshotIdentifier
          .replace("%s", queryId)
          .replace("%o", organismId)}`,
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

  const queryUrlSuffix = urlSuffix
    .replace("%s", queryId)
    .replace("%o", organismId);
  const queryUrl = queryLatestSnapshotIdentifier
    ? `${CELLGUIDE_DATA_URL_WITH_RDEV_SUFFIX}/${latestSnapshotIdentifier}/${queryUrlSuffix}`
    : `${CELLGUIDE_DATA_URL_WITH_RDEV_SUFFIX}/${queryUrlSuffix}`;

  return useQuery(
    queryId || organismId
      ? [queryKey, queryId, organismId, latestSnapshotIdentifier]
      : [queryKey],
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

/* ========== valid_explorer_cxgs ========== */
export const USE_VALID_EXPLORER_CXGS_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_VALID_EXPLORER_CXGS],
  id: "cell-guide-valid-explorer-cxgs-query",
};

export interface ValidExplorerCxgsQueryResponse {
  organism_celltype_cxgs: { [organism: string]: string[] };
  organism_tissue_celltype_cxgs: {
    [organism: string]: { [tissue: string]: string[] };
  };
}

export const useValidExplorerCxgs = () => {
  return useCellGuideQuery<ValidExplorerCxgsQueryResponse>({
    dataType: TYPES.VALID_EXPLORER_CXGS,
  });
};
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

export const useCellOntologyTree = (
  organismName: string
): UseQueryResult<CellOntologyTreeResponse> => {
  return useCellGuideQuery<CellOntologyTreeResponse>({
    dataType: TYPES.CELL_ONTOLOGY_TREE,
    organismName,
  });
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
  entityId: string,
  organismName: string
): UseQueryResult<CellOntologyTreeStateResponse> => {
  return useCellGuideQuery<CellOntologyTreeStateResponse>({
    dataType: TYPES.CELL_ONTOLOGY_TREE_STATE_CELLTYPE,
    queryId: entityId,
    organismName,
  });
};

/* ========== ontology_tree_state_tissue ========== */
export const USE_CELL_ONTOLOGY_TREE_STATE_TISSUE_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_CELL_ONTOLOGY_TREE_STATE_TISSUE],
  id: "cell-guide-cell-ontology-tree-state-tissue-query",
};

export const useCellOntologyTreeStateTissue = (
  entityId: string,
  organismName: string
): UseQueryResult<CellOntologyTreeStateResponse> => {
  return useCellGuideQuery<CellOntologyTreeStateResponse>({
    dataType: TYPES.CELL_ONTOLOGY_TREE_STATE_TISSUE,
    queryId: entityId,
    organismName,
  });
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
  return useCellGuideQuery<SourceCollectionsQueryResponse>({
    dataType: TYPES.SOURCE_COLLECTIONS,
    queryId: entityId,
  });
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
  specificity: number;
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
  return useCellGuideQuery<ComputationalMarkersQueryResponse>({
    dataType: TYPES.COMPUTATIONAL_MARKERS,
    queryId: entityId,
  });
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
  return useCellGuideQuery<CanonicalMarkersQueryResponse>({
    dataType: TYPES.CANONICAL_MARKERS,
    queryId: entityId,
  });
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
  return useCellGuideQuery<string>({
    dataType: TYPES.GPT_DESCRIPTION,
    queryId: entityId,
    queryLatestSnapshotIdentifier: false,
  });
};

/* ========== description ========== */
export const USE_VALIDATED_DESCRIPTION_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_VALIDATED_DESCRIPTION],
  id: "cell-guide-validated-description-query",
};

interface ValidatedDescriptionQueryResponse {
  description: string;
  references: string[];
}

export const useValidatedDescription = (
  entityId: string
): UseQueryResult<ValidatedDescriptionQueryResponse> => {
  return useCellGuideQuery<ValidatedDescriptionQueryResponse>({
    dataType: TYPES.VALIDATED_DESCRIPTION,
    queryId: entityId,
    queryLatestSnapshotIdentifier: false,
  });
};

/* ========== SEO description ========== */
export const USE_GPT_SEO_DESCRIPTION_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_GPT_SEO_DESCRIPTION],
  id: "cell-guide-gpt-seo-description-query",
};

export const fetchGptSeoDescription = async (
  entityId: string
): Promise<string> => {
  entityId = entityId.replace(":", "_");
  // This function is used server-side to fetch the GPT SEO description.
  const url = `${CELLGUIDE_DATA_URL_WITH_RDEV_SUFFIX}/${QUERY_MAPPING[
    TYPES.GPT_SEO_DESCRIPTION
  ].urlSuffix.replace("%s", entityId)}`;
  const response = await fetch(url);
  if (!response.ok) {
    throw new Error(`HTTP error! status: ${response.status}`);
  }
  return await response.json();
};

/* ========== marker_gene_presence ========== */
export const USE_MARKER_GENE_PRESENCE_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_MARKER_GENE_PRESENCE],
  id: "cell-guide-marker-gene-presence-query",
};

interface MarkerGenePresenceQueryResponse {
  [gene: string]: {
    [organism: string]: {
      [tissue: string]: {
        me: number;
        pc: number;
        marker_score: number;
        specificity: number;
        cell_type_id: string;
      }[];
    };
  };
}

export const useMarkerGenePresenceQuery =
  (): UseQueryResult<MarkerGenePresenceQueryResponse> => {
    return useCellGuideQuery<MarkerGenePresenceQueryResponse>({
      dataType: TYPES.MARKER_GENE_PRESENCE,
    });
  };

export const USE_CELLTYPE_TISSUE_MAPPING_QUERY = {
  entities: [ENTITIES.CELL_GUIDE_CELLTYPE_TISSUE_MAPPING],
  id: "cell-guide-celltype-tissue-mapping-query",
};

interface CellTypeTissueMappingResponse {
  [cellTypeId: string]: string[];
}

export const useCellTypeTissueMapping = (
  organismName: string
): UseQueryResult<CellTypeTissueMappingResponse> => {
  return useCellGuideQuery<CellTypeTissueMappingResponse>({
    dataType: TYPES.CELLTYPE_TISSUE_MAPPING,
    organismName,
  });
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
    return useCellGuideQuery<CellTypeMetadataQueryResponse>({
      dataType: TYPES.CELLTYPE_METADATA,
    });
  };

// this naked fetch is used for SSR
export const fetchCellTypeMetadata =
  async (): Promise<CellTypeMetadataQueryResponse> => {
    // This function is used server-side to fetch the GPT SEO description.
    const latestSnapshotIdentifierUrl = `${CELLGUIDE_DATA_URL_WITH_RDEV_SUFFIX}/${
      QUERY_MAPPING[TYPES.LATEST_SNAPSHOT_IDENTIFIER].urlSuffix
    }`;
    const latestSnapshotIdentifierResponse = await fetch(
      latestSnapshotIdentifierUrl
    );
    const latestSnapshotIdentifier =
      await latestSnapshotIdentifierResponse.text();
    const url = `${CELLGUIDE_DATA_URL_WITH_RDEV_SUFFIX}/${latestSnapshotIdentifier}/${
      QUERY_MAPPING[TYPES.CELLTYPE_METADATA].urlSuffix
    }`;
    const response = await fetch(url);
    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }
    return await response.json();
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
    return useCellGuideQuery<TissueMetadataQueryResponse>({
      dataType: TYPES.TISSUE_METADATA,
    });
  };

export const fetchTissueMetadata =
  async (): Promise<TissueMetadataQueryResponse> => {
    // This function is used server-side to fetch the GPT SEO description.
    const latestSnapshotIdentifierUrl = `${CELLGUIDE_DATA_URL_WITH_RDEV_SUFFIX}/${
      QUERY_MAPPING[TYPES.LATEST_SNAPSHOT_IDENTIFIER].urlSuffix
    }`;
    const latestSnapshotIdentifierResponse = await fetch(
      latestSnapshotIdentifierUrl
    );
    const latestSnapshotIdentifier =
      await latestSnapshotIdentifierResponse.text();
    const url = `${CELLGUIDE_DATA_URL_WITH_RDEV_SUFFIX}/${latestSnapshotIdentifier}/${
      QUERY_MAPPING[TYPES.TISSUE_METADATA].urlSuffix
    }`;
    const response = await fetch(url);
    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }
    return await response.json();
  };

/* ========== Lookup tables for tissues ========== */
export function useAllTissuesLookupTables(cellTypeId: string): {
  allTissuesLabelToIdLookup: Map<string, string>;
  computationalMarkers: ComputationalMarkersQueryResponse;
} {
  const { data: tissueData } = useTissueMetadata();
  const { data: computationalMarkers } = useComputationalMarkers(cellTypeId);

  return useMemo(() => {
    if (!tissueData || !computationalMarkers) {
      return {
        allTissuesLabelToIdLookup: new Map<string, string>(),
        computationalMarkers: EMPTY_ARRAY,
      };
    }
    const tissueIdByLabel: { [label: string]: string } = {};

    for (const tissueId in tissueData) {
      const tissue = tissueData[tissueId];
      tissueIdByLabel[tissue.name] = tissue.id;
    }

    const allTissuesLabelToIdLookup = new Map<string, string>();

    for (const markerGene of computationalMarkers) {
      const label = markerGene.groupby_dims.tissue_ontology_term_label;
      if (!label) continue;
      allTissuesLabelToIdLookup.set(label, tissueIdByLabel[label]);
    }
    return { allTissuesLabelToIdLookup, computationalMarkers };
  }, [tissueData, computationalMarkers]);
}

/**
 * Mapping from data/response type to properties used for querying
 */
const QUERY_MAPPING: {
  [key in TYPES]: CellGuideQuery;
} = {
  VALID_EXPLORER_CXGS: {
    queryKey: USE_VALID_EXPLORER_CXGS_QUERY,
    urlSuffix: "valid_explorer_cxgs.json",
  },
  CELL_ONTOLOGY_TREE: {
    queryKey: USE_CELL_ONTOLOGY_TREE_QUERY,
    urlSuffix: "ontology_tree/%o/ontology_graph.json",
  },
  CELL_ONTOLOGY_TREE_STATE_CELLTYPE: {
    queryKey: USE_CELL_ONTOLOGY_TREE_STATE_CELLTYPE_QUERY,
    urlSuffix: `ontology_tree/%o/cell_type_ontology_tree_state/%s.json`,
  },
  CELL_ONTOLOGY_TREE_STATE_TISSUE: {
    queryKey: USE_CELL_ONTOLOGY_TREE_STATE_TISSUE_QUERY,
    urlSuffix: `ontology_tree/%o/tissue_ontology_tree_state/%s.json`,
  },
  SOURCE_COLLECTIONS: {
    queryKey: USE_SOURCE_COLLECTIONS_QUERY,
    urlSuffix: `source_collections/%s.json`,
  },
  COMPUTATIONAL_MARKERS: {
    queryKey: USE_COMPUTATIONAL_MARKERS_QUERY,
    urlSuffix: `computational_marker_genes/%s.json`,
  },
  MARKER_GENE_PRESENCE: {
    queryKey: USE_MARKER_GENE_PRESENCE_QUERY,
    urlSuffix: `computational_marker_genes/marker_gene_presence.json.gz`,
  },
  CANONICAL_MARKERS: {
    queryKey: USE_CANONICAL_MARKERS_QUERY,
    urlSuffix: `canonical_marker_genes/%s.json`,
  },
  GPT_DESCRIPTION: {
    queryKey: USE_GPT_DESCRIPTION_QUERY,
    urlSuffix: `gpt_descriptions/%s.json`,
  },
  VALIDATED_DESCRIPTION: {
    queryKey: USE_VALIDATED_DESCRIPTION_QUERY,
    urlSuffix: `validated_descriptions/%s.json`,
  },
  GPT_SEO_DESCRIPTION: {
    queryKey: USE_GPT_SEO_DESCRIPTION_QUERY,
    urlSuffix: `gpt_seo_descriptions/%s.json`,
  },
  CELLTYPE_METADATA: {
    queryKey: USE_CELLTYPE_METADATA_QUERY,
    urlSuffix: "celltype_metadata.json",
  },
  CELLTYPE_TISSUE_MAPPING: {
    queryKey: USE_CELLTYPE_TISSUE_MAPPING_QUERY,
    urlSuffix: "ontology_tree/%o/celltype_to_tissue_mapping.json",
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
