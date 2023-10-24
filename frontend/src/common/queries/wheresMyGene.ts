import { useContext, useMemo } from "react";
import { useQuery, UseQueryResult } from "react-query";
import { API_URL } from "src/configs/configs";
import { FMG_GENE_STRENGTH_THRESHOLD } from "src/views/WheresMyGeneV2/common/constants";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGeneV2/common/store";
import { setSnapshotId } from "src/views/WheresMyGeneV2/common/store/actions";
import {
  CellType,
  CellTypeGeneExpressionSummaryData,
  CellTypeId,
  CellTypeSummary,
  CompareOptionId,
  GeneExpressionSummary,
  GeneInfo,
  RawCellTypeGeneExpressionSummaryData,
  ViewId,
  Organism as IOrganism,
} from "src/views/WheresMyGeneV2/common/types";
import { API } from "../API";
import { APIV2 } from "src/common/tempAPIV2";

import { ROUTES } from "../constants/routes";
import { EMPTY_OBJECT } from "../constants/utils";
import { DEFAULT_FETCH_OPTIONS, JSON_BODY_FETCH_OPTIONS } from "./common";
import { ENTITIES } from "./entities";
import { Dataset } from "@mui/icons-material";
import { formatCitation } from "../utils/formatCitation";

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

export interface OntologyTerm {
  id: string;
  name: string;
}
interface OntologyTermsByOrganism {
  [organismID: string]: Array<OntologyTerm>;
}

export interface PrimaryFilterDimensionsResponse {
  genes: OntologyTermsByOrganism;
  organisms: Array<OntologyTerm>;
  snapshotId: string;
  tissues: OntologyTermsByOrganism;
}

function replaceV1WithV2(version: 1 | 2) {
  return version === 1 ? API : APIV2;
}

export async function fetchPrimaryFilterDimensions(
  version: 1 | 2
): Promise<PrimaryFilterDimensionsResponse> {
  const url = API_URL + replaceV1WithV2(version).WMG_PRIMARY_FILTER_DIMENSIONS;

  const response: RawPrimaryFilterDimensionsResponse = await (
    await fetch(url, DEFAULT_FETCH_OPTIONS)
  ).json();

  return transformPrimaryFilterDimensions(response);
}

function flattenOntologyTermsByOrganism(
  termsObject: RawOntologyTermsByOrganism
): OntologyTermsByOrganism {
  return Object.entries(termsObject).reduce((memo, [organismId, genes]) => {
    memo[organismId] = genes.map(toEntity);
    return memo;
  }, {} as OntologyTermsByOrganism);
}

export function generateTermsByKey(
  flattenedTerms: OntologyTermsByOrganism,
  key: keyof OntologyTerm
): {
  [key: string]: OntologyTerm;
} {
  const termsByKey: { [key: string]: OntologyTerm } = {};

  Object.values(flattenedTerms).forEach((terms) => {
    for (const term of terms) {
      termsByKey[term[key]] = term;
    }
  });

  return termsByKey;
}

function transformPrimaryFilterDimensions(
  response: RawPrimaryFilterDimensionsResponse
): PrimaryFilterDimensionsResponse {
  const { gene_terms, organism_terms, snapshot_id, tissue_terms } = response;

  return {
    genes: flattenOntologyTermsByOrganism(gene_terms),
    organisms: organism_terms.map(toEntity),
    snapshotId: snapshot_id,
    tissues: flattenOntologyTermsByOrganism(tissue_terms),
  };
}

export const USE_PRIMARY_FILTER_DIMENSIONS = {
  entities: [ENTITIES.WMG_PRIMARY_FILTER_DIMENSIONS],
  id: "wmg-primaryFilterDimensions",
};

export function usePrimaryFilterDimensions(
  version: 1 | 2 = 2
): UseQueryResult<PrimaryFilterDimensionsResponse> {
  const dispatch = useContext(DispatchContext);

  // (thuang): Refresh query when the snapshotId changes
  const currentSnapshotId = useSnapshotId();

  return useQuery<PrimaryFilterDimensionsResponse>(
    [USE_PRIMARY_FILTER_DIMENSIONS, currentSnapshotId, version],
    () => fetchPrimaryFilterDimensions(version),
    {
      onSuccess(response) {
        if (!response || !dispatch) return;

        const { snapshotId } = response;

        if (currentSnapshotId !== snapshotId) {
          dispatch(setSnapshotId(snapshotId));
        }
      },
      // (thuang): We don't need to refetch during the session
      staleTime: Infinity,
    }
  );
}

const TEMP_ALLOW_NAME_LIST = ["Homo sapiens", "Mus musculus"];

export function useAvailableOrganisms(version: 1 | 2 = 2) {
  const { data, isLoading } = usePrimaryFilterDimensions(version);

  if (isLoading) {
    return { isLoading, data: null };
  }

  return {
    isLoading,
    data: data?.organisms.filter((organism: IOrganism) =>
      TEMP_ALLOW_NAME_LIST.includes(organism.name)
    ),
  };
}

interface FilterV1 {
  gene_ontology_term_ids: string[];
  organism_ontology_term_id: string;
  tissue_ontology_term_ids: string[];
  dataset_ids: string[];
  disease_ontology_term_ids: string[];
  sex_ontology_term_ids: string[];
  development_stage_ontology_term_ids: string[];
  self_reported_ethnicity_ontology_term_ids: string[];
}

interface FilterV2 {
  gene_ontology_term_ids: string[];
  organism_ontology_term_id: string;
  dataset_ids: string[];
  disease_ontology_term_ids: string[];
  sex_ontology_term_ids: string[];
  development_stage_ontology_term_ids: string[];
  self_reported_ethnicity_ontology_term_ids: string[];
}

type Filter = FilterV1 | FilterV2;

interface FilterSecondary {
  organism_ontology_term_id: string;
  tissue_ontology_term_ids?: string[];
  dataset_ids: string[];
  disease_ontology_term_ids: string[];
  sex_ontology_term_ids: string[];
  development_stage_ontology_term_ids: string[];
  self_reported_ethnicity_ontology_term_ids: string[];
}

export interface FiltersQuery {
  filter: FilterSecondary;
}

export interface Query {
  filter: Filter;
  compare?: CompareOptionId;
}

export type OptionId = string;

export interface QueryResponse {
  expression_summary: {
    // gene_ontology_term_id
    [geneId: string]: {
      [tissueId: string]: {
        tissue_stats: {
          aggregated: RawCellTypeGeneExpressionSummaryData;
        };
      } & {
        [cellTypeId: CellTypeId]: {
          aggregated: RawCellTypeGeneExpressionSummaryData;
          [
            compareOptionId: CompareOptionId
          ]: RawCellTypeGeneExpressionSummaryData;
        };
      };
    };
  };
  snapshot_id: string;
  term_id_labels: {
    cell_types: {
      [tissue_type_ontology_term_id: string]: TissueStats & CellTypeStats;
    };
    genes: {
      [id: string]: string;
    }[];
  };
}

interface TissueStats {
  tissue_stats: {
    aggregated: {
      name: string;
      order: -1;
      tissue_ontology_term_id: CellTypeId;
      total_count: number;
    };
  };
}

interface CellTypeStats {
  [cell_type_ontology_term_id: string]: {
    [compareOptionId: CompareOptionId]: {
      name: string;
      cell_type_ontology_term_id: CellTypeId;
      order: number;
      total_count: number;
    };
  };
}

function isTissueStats(
  stats:
    | TissueStats["tissue_stats"]["aggregated"]
    | CellTypeStats[string][CompareOptionId]
): stats is TissueStats["tissue_stats"]["aggregated"] {
  return (
    (stats as TissueStats["tissue_stats"]["aggregated"])
      .tissue_ontology_term_id !== undefined
  );
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
    publication_citations: string[];
    development_stage_terms: { [id: string]: string }[];
    self_reported_ethnicity_terms: { [id: string]: string }[];
    tissue_terms: { [id: string]: string }[];
  };
  snapshot_id: string;
}

async function fetchFiltersQuery({
  query,
  signal,
  version,
}: {
  query: FiltersQuery | null;
  signal?: AbortSignal;
  version: 1 | 2;
}): Promise<FiltersQueryResponse | undefined> {
  if (!query) return;

  const url = API_URL + replaceV1WithV2(version).WMG_FILTERS_QUERY;

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

async function fetchQuery({
  query,
  signal,
  version,
}: {
  query: Query | null;
  signal?: AbortSignal;
  version: 1 | 2;
}): Promise<QueryResponse | undefined> {
  if (!query) return;

  const url = API_URL + replaceV1WithV2(version).WMG_QUERY;

  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    body: JSON.stringify(query),
    method: "POST",
    signal,
  });
  const json: QueryResponse = await response.json();

  if (!response.ok) {
    throw json;
  }

  return json;
}

export const USE_QUERY = {
  entities: [ENTITIES.WMG_QUERY],
  id: "wmg-query",
};

export const USE_FILTERS_QUERY = {
  entities: [ENTITIES.WMG_FILTERS_QUERY],
  id: "wmg-filters-query",
};

export function useWMGQuery(
  query: Query | null,
  version: 1 | 2 = 2
): UseQueryResult<QueryResponse> {
  const dispatch = useContext(DispatchContext);

  // (thuang): Refresh query when the snapshotId changes
  const currentSnapshotId = useSnapshotId();

  query = clobberQueryIfSubsetOfPrev(query, [
    "gene_ontology_term_ids",
    "tissue_ontology_term_ids",
  ]);

  return useQuery(
    [USE_QUERY, query, currentSnapshotId],
    ({ signal }) => fetchQuery({ query, signal, version }),
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

export function useWMGFiltersQuery(
  query: FiltersQuery | null,
  version: 1 | 2 = 2
): UseQueryResult<FiltersQueryResponse> {
  const dispatch = useContext(DispatchContext);

  // (thuang): Refresh query when the snapshotId changes
  const currentSnapshotId = useSnapshotId();

  return useQuery(
    [USE_FILTERS_QUERY, query, currentSnapshotId],
    ({ signal }) => fetchFiltersQuery({ query, signal, version }),
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
  publication_citations: [],
  sex_terms: [],
  tissue_terms: [],
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
  publication_citations: { id: string; name: string }[];
  sex_terms: { id: string; name: string }[];
  tissue_terms: { id: string; name: string }[];
}

export function useFilterDimensions(version: 1 | 2 = 2): {
  data: FilterDimensions;
  isLoading: boolean;
} {
  const requestBody = useWMGFiltersQueryRequestBody(version);
  const { data, isLoading } = useWMGFiltersQuery(requestBody);

  return useMemo(() => {
    if (isLoading || !data) return { data: EMPTY_FILTER_DIMENSIONS, isLoading };

    const { filter_dims } = data;

    const {
      datasets,
      development_stage_terms,
      disease_terms,
      self_reported_ethnicity_terms,
      publication_citations,
      sex_terms,
      tissue_terms,
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
        publication_citations: publication_citations.map(toCitationEntity),
        sex_terms: sex_terms.map(toEntity),
        tissue_terms: tissue_terms.map(toEntity),
      },
      isLoading: false,
    };
  }, [data, isLoading]);
}

export function useExpressionSummary(version: 1 | 2 = 2): {
  isLoading: boolean;
  data: QueryResponse["expression_summary"];
} {
  const requestBody = useWMGQueryRequestBody(version);
  const { data, isLoading } = useWMGQuery(requestBody, version);

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
  [tissueName: string]: CellTypeRow[];
}

// Cell types that we want to exclude from each tissue

const FILTERED_CELL_TYPE_ONTOLOGY_IDS = [
  "CL:0000003", // Native cell
  "CL:0000255", // Eukaryotic cell
  "CL:0000548", // Animal cell
];

export function useCellTypesByTissueName(version: 1 | 2 = 2): {
  isLoading: boolean;
  data: CellTypeByTissueName;
} {
  const { isLoading } = useExpressionSummary(version);

  const {
    data: primaryFilterDimensions,
    isLoading: isLoadingPrimaryFilterDimensions,
  } = usePrimaryFilterDimensions(version);

  const { data: termIdLabels, isLoading: isLoadingTermIdLabels } =
    useTermIdLabels(version);

  return useMemo(() => {
    if (
      isLoading ||
      isLoadingPrimaryFilterDimensions ||
      !primaryFilterDimensions ||
      isLoadingTermIdLabels ||
      !Object.keys(termIdLabels.cell_types).length
    ) {
      return { data: EMPTY_OBJECT, isLoading };
    }

    const { tissues } = primaryFilterDimensions;

    const tissuesById = generateTermsByKey(tissues, "id");

    const cellTypesByTissueName: { [tissueName: string]: CellType[] } = {};

    for (const [tissueId, rawTissueCellTypes] of Object.entries(
      termIdLabels.cell_types
    )) {
      const tissueName = tissuesById[tissueId].name;

      const tissueCellTypes = Object.entries(rawTissueCellTypes)
        // (thuang): Reverse the order, so the first cell type is at the top of
        // the heat map
        .reverse()
        .filter(([_, cellTypeRow]) => {
          return !FILTERED_CELL_TYPE_ONTOLOGY_IDS.includes(cellTypeRow.id);
        })
        .map(([_, cellTypeName]) => {
          return cellTypeName;
        });

      cellTypesByTissueName[tissueName] = tissueCellTypes;
    }

    return {
      data: cellTypesByTissueName,
      isLoading,
    };
  }, [
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

export function useGeneExpressionSummariesByTissueName(version: 1 | 2 = 2): {
  data: GeneExpressionSummariesByTissueName;
  isLoading: boolean;
} {
  const { data, isLoading } = useExpressionSummary(version);
  const {
    data: primaryFilterDimensions,
    isLoading: isLoadingPrimaryFilterDimensions,
  } = usePrimaryFilterDimensions(version);

  const { data: termIdLabels, isLoading: isLoadingTermIdLabels } =
    useTermIdLabels(version);

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
    const tissuesById = generateTermsByKey(tissues, "id");
    const result: GeneExpressionSummariesByTissueName = {};

    for (const [geneId, expressionSummariesByTissue] of Object.entries(data)) {
      const geneName = termIdLabels.genes[geneId];

      for (const [tissueId, expressionSummariesByCellType] of Object.entries(
        expressionSummariesByTissue
      )) {
        const mergedExpressionSummaries = mergeExpressionSummaries(
          expressionSummariesByCellType,
          tissueId
        );

        updateResult(
          result,
          tissuesById[tissueId].name,
          geneName,
          mergedExpressionSummaries
        );
      }
    }

    return { data: result, isLoading };

    /**
     * This produces a collection that contains all the expression summaries for
     * a given tissue, including all the cell types and compare options
     */
    function mergeExpressionSummaries(
      expressionSummariesByCellType: QueryResponse["expression_summary"][string][string],
      tissueId: string
    ) {
      const mergedExpressionSummaries = [];

      for (const [
        cellTypeId,
        expressionSummariesByCompareOption,
      ] of Object.entries(expressionSummariesByCellType)) {
        for (const [compareOptionId, expressionSummary] of Object.entries(
          expressionSummariesByCompareOption
        )) {
          mergedExpressionSummaries.push(
            transformExpressionSummary(
              expressionSummary,
              cellTypeId,
              tissueId,
              compareOptionId
            )
          );
        }
      }

      return mergedExpressionSummaries;
    }

    function transformExpressionSummary(
      expressionSummary: RawCellTypeGeneExpressionSummaryData,
      cellTypeId: string,
      tissueId: string,
      compareOptionId: string
    ) {
      return transformCellTypeGeneExpressionSummaryData({
        ...expressionSummary,
        viewId: getCellTypeViewId(
          cellTypeId === "tissue_stats" ? tissueId : cellTypeId,
          compareOptionId
        ),
      });
    }

    function updateResult(
      result: GeneExpressionSummariesByTissueName,
      tissueName: string,
      geneName: string,
      mergedExpressionSummaries: CellTypeGeneExpressionSummaryData[]
    ) {
      const tissueGeneExpressionSummaries = result[tissueName] || {};

      tissueGeneExpressionSummaries[geneName] = {
        cellTypeGeneExpressionSummaries: mergedExpressionSummaries,
        name: geneName,
      };

      result[tissueName] = tissueGeneExpressionSummaries;
    }
  }, [
    data,
    isLoading,
    primaryFilterDimensions,
    isLoadingPrimaryFilterDimensions,
    termIdLabels,
    isLoadingTermIdLabels,
  ]);
}

type TransformCellTypeGeneExpressionSummaryDataInput =
  RawCellTypeGeneExpressionSummaryData & {
    viewId: CellTypeGeneExpressionSummaryData["viewId"];
  };

function transformCellTypeGeneExpressionSummaryData(
  data: TransformCellTypeGeneExpressionSummaryDataInput
): CellTypeGeneExpressionSummaryData {
  const { viewId, pc, me, tpc, n } = data;

  return {
    ...data,
    expressedCellCount: n,
    meanExpression: me,
    percentage: pc,
    tissuePercentage: tpc,
    viewId,
  };
}

interface TermIdLabels {
  cell_types: {
    [tissueID: string]: {
      [viewId: ViewId]: {
        id: string;
        name: string;
        order: number;
        total_count: number;
        viewId: ViewId;
        isAggregated: boolean;
      };
    };
  };
  genes: { [id: string]: string };
}

/**
 * (thuang): If BE changes the aggregated option id from "aggregated" to
 * something else, we'll need to update it here
 */
export const COMPARE_OPTION_ID_FOR_AGGREGATED = "aggregated";

/**
 * (thuang): If BE changes the unknown option id from "unknown" to
 * something else, we'll need to update it here
 */
export const COMPARE_OPTION_ID_FOR_UNKNOWN = "unknown";

export interface CellTypeRow {
  name: CellTypeSummary["name"];
  order: CellTypeSummary["order"];
  total_count: CellTypeSummary["total_count"];
  viewId: CellTypeSummary["viewId"];
  id: CellTypeSummary["id"];
  // (thuang): boolean flag if the cell type is an aggregated option or not
  isAggregated: boolean;
  cellTypeName: string;
}

export function useTermIdLabels(version: 1 | 2 = 2): {
  data: TermIdLabels;
  isLoading: boolean;
} {
  const requestBody = useWMGQueryRequestBody(version);

  const { data, isLoading } = useWMGQuery(requestBody, version);

  return useMemo(() => {
    if (isLoading || !data) {
      return {
        data: { cell_types: EMPTY_OBJECT, genes: EMPTY_OBJECT },
        isLoading,
      };
    }

    const {
      term_id_labels: { cell_types, genes },
    } = data;

    const returnCellTypes: TermIdLabels["cell_types"] = {};

    Object.entries(cell_types).forEach(
      ([tissueID, tissueCellTypesWithCompareOptions]) => {
        const sortedTissueCellTypesWithCompareOptions =
          getSortedTissueCellTypesWithCompareOptions(
            tissueCellTypesWithCompareOptions
          );

        const result: {
          [viewId: CellTypeRow["viewId"]]: CellTypeRow;
        } = {};
        if (tissueCellTypesWithCompareOptions["tissue_stats"])
          addCellTypeRowToResult({
            result,
            sortedCellTypeCompareOptions: [
              [
                "aggregated",
                tissueCellTypesWithCompareOptions["tissue_stats"]["aggregated"],
              ],
            ],
          });

        for (const cellTypeWithCompareOptions of sortedTissueCellTypesWithCompareOptions) {
          const sortedCellTypeCompareOptions = getSortedCellTypeCompareOptions(
            cellTypeWithCompareOptions
          );

          addCellTypeRowToResult({
            result,
            sortedCellTypeCompareOptions,
          });
        }

        returnCellTypes[tissueID] = result;
      }
    );

    return {
      data: {
        cell_types: returnCellTypes,
        genes: aggregateIdLabels(genes),
      },
      isLoading: false,
    };
  }, [data, isLoading]);
}

/**
 * (thuang): This function sorts a tissue's cellTypeWithCompareOptions objects
 * by `order`
 */
function getSortedTissueCellTypesWithCompareOptions({
  tissue_stats: _,
  ...tissueCellTypesWithCompareOptions
}: QueryResponse["term_id_labels"]["cell_types"][string]): QueryResponse["term_id_labels"]["cell_types"][string][string][] {
  return Object.values(tissueCellTypesWithCompareOptions).sort((a, b) => {
    const aAggregated = a.aggregated;
    const bAggregated = b.aggregated;

    const aOrder = aAggregated.order;
    const bOrder = bAggregated.order;

    if (aOrder < bOrder) {
      return -1;
    }

    if (aOrder > bOrder) {
      return 1;
    }

    return 0;
  });
}

/**
 * (thuang): This sorts all the compare options for a cellType by the order:
 * aggregated, compareOptionId, unknown
 */
function getSortedCellTypeCompareOptions(
  cellTypeWithCompareOptions: QueryResponse["term_id_labels"]["cell_types"][string][string]
): [
  // compareOptionId
  string,
  QueryResponse["term_id_labels"]["cell_types"][string][string][string],
][] {
  return Object.entries(cellTypeWithCompareOptions).sort((a, b) => {
    const aCompareOptionId = a[0];
    const bCompareOptionId = b[0];

    const aIsAggregated = aCompareOptionId === COMPARE_OPTION_ID_FOR_AGGREGATED;
    const bIsAggregated = bCompareOptionId === COMPARE_OPTION_ID_FOR_AGGREGATED;

    const aIsUnknown = aCompareOptionId === COMPARE_OPTION_ID_FOR_UNKNOWN;
    const bIsUnknown = bCompareOptionId === COMPARE_OPTION_ID_FOR_UNKNOWN;

    if (aIsAggregated) {
      return -1;
    }

    if (bIsAggregated) {
      return 1;
    }

    if (aIsUnknown) {
      return 1;
    }

    if (bIsUnknown) {
      return -1;
    }

    const aCompareOption = a[1];
    const bCompareOption = b[1];

    // cchoi: higher number goes first!
    return bCompareOption.total_count - aCompareOption.total_count;
  });
}

function addCellTypeRowToResult({
  result,
  sortedCellTypeCompareOptions,
}: {
  result: {
    [viewId: CellTypeRow["viewId"]]: CellTypeRow;
  };
  sortedCellTypeCompareOptions: [
    // compareOptionId
    string,
    (
      | CellTypeStats[string][CompareOptionId]
      | TissueStats["tissue_stats"]["aggregated"]
    ),
  ][];
}) {
  let cellTypeName = "";

  for (const [
    compareOptionId,
    compareOptionData,
  ] of sortedCellTypeCompareOptions) {
    const isAggregated = compareOptionId === COMPARE_OPTION_ID_FOR_AGGREGATED;

    const { name: rawName, total_count, order } = compareOptionData;

    const termIsTissueStats = isTissueStats(compareOptionData);

    const termID = isTissueStats(compareOptionData)
      ? compareOptionData.tissue_ontology_term_id
      : compareOptionData.cell_type_ontology_term_id;

    /**
     * (thuang): We manually indent 4 spaces instead of using CSS, so we don't
     * have to update SVG render function to mimic CSS padding
     */
    const name = isAggregated
      ? rawName
      : compareOptionId === COMPARE_OPTION_ID_FOR_UNKNOWN
      ? `    ${rawName || COMPARE_OPTION_ID_FOR_UNKNOWN}`
      : `    ${rawName || compareOptionId}`;

    if (isAggregated) {
      cellTypeName = rawName;
    }

    const viewId = getCellTypeViewId(termID, compareOptionId);

    result[viewId] = {
      cellTypeName: termIsTissueStats ? termID : cellTypeName,
      id: termID,
      isAggregated,
      name,
      order,
      total_count,
      viewId,
    };
  }
}

function aggregateIdLabels(items: { [id: string]: string }[]): {
  [id: string]: string;
} {
  return items.reduce((memo, item) => ({ ...memo, ...item }), {});
}

function useWMGQueryRequestBody(version: 1 | 2) {
  const {
    compare,
    selectedGenes,
    selectedTissues,
    selectedOrganismId,
    selectedFilters,
  } = useContext(StateContext);

  const { data } = usePrimaryFilterDimensions(version);

  const {
    datasets,
    developmentStages,
    diseases,
    ethnicities,
    sexes,
    publications,
  } = selectedFilters;

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

  // seve: This is a nice mapping that may start being used in a few places across WMG, might be worth moving to a global store.
  const tissuesByName = useMemo(() => {
    let result: { [name: string]: OntologyTerm } = {};

    if (!data) return result;

    const { tissues } = data;

    result = generateTermsByKey(tissues, "name");

    return result;
  }, [data]);

  return useMemo(() => {
    if (
      !data ||
      !selectedOrganismId ||
      (version === 1 && !selectedTissues?.length)
    ) {
      return null;
    }
    const gene_ontology_term_ids = selectedGenes.map((geneName) => {
      return organismGenesByName[geneName].id;
    });
    if (!gene_ontology_term_ids.length) gene_ontology_term_ids.push(".");
    const tissue_ontology_term_ids = selectedTissues?.map((tissueName) => {
      return tissuesByName[tissueName].id;
    });

    return {
      compare,
      filter: {
        dataset_ids: datasets,
        development_stage_ontology_term_ids: developmentStages,
        disease_ontology_term_ids: diseases,
        gene_ontology_term_ids,
        organism_ontology_term_id: selectedOrganismId,
        self_reported_ethnicity_ontology_term_ids: ethnicities,
        sex_ontology_term_ids: sexes,
        publication_citations: publications,
        ...(version === 1 && { tissue_ontology_term_ids }),
      },
      is_rollup: true, // this could be made toggleable by users in the future
    };
  }, [
    data,
    selectedOrganismId,
    version,
    selectedTissues,
    selectedGenes,
    compare,
    datasets,
    developmentStages,
    diseases,
    ethnicities,
    sexes,
    organismGenesByName,
    tissuesByName,
    publications,
  ]);
}

function useWMGFiltersQueryRequestBody(
  version: 1 | 2 = 2
): FiltersQuery | null {
  const { selectedOrganismId, selectedFilters, filteredCellTypeIds } =
    useContext(StateContext);

  const { data } = usePrimaryFilterDimensions(version);

  const {
    tissues,
    datasets,
    developmentStages,
    diseases,
    ethnicities,
    sexes,
    publications,
  } = selectedFilters;

  return useMemo(() => {
    if (!data || !selectedOrganismId) {
      return null;
    }

    return {
      filter: {
        dataset_ids: datasets,
        development_stage_ontology_term_ids: developmentStages,
        disease_ontology_term_ids: diseases,
        organism_ontology_term_id: selectedOrganismId,
        self_reported_ethnicity_ontology_term_ids: ethnicities,
        sex_ontology_term_ids: sexes,
        tissue_ontology_term_ids: tissues,
        publication_citations: publications,
        cell_type_ontology_term_ids: filteredCellTypeIds,
      },
    };
  }, [
    tissues,
    selectedOrganismId,
    data,
    datasets,
    developmentStages,
    diseases,
    ethnicities,
    publications,
    sexes,
    filteredCellTypeIds,
  ]);
}

let prevQuery: Query | null;

type KeysOfUnion<T> = T extends T ? keyof T : never;

function clobberQueryIfSubsetOfPrev(
  query: Query | null,
  filtersToCheck: KeysOfUnion<Filter>[]
): Query | null {
  if (prevQuery == query) return prevQuery;
  if (!prevQuery || !query) {
    prevQuery = query;
    return query;
  }

  if (prevQuery.compare !== query.compare) return query;

  if (
    (Object.entries(query.filter) as [keyof Filter, string[]][]).every(
      ([key, value]) => {
        //just check for equality on the filters we aren't checking for subsets
        if (!filtersToCheck.includes(key))
          return (
            JSON.stringify(value) === JSON.stringify(prevQuery?.filter[key])
          );
        return value.every((elem) => prevQuery?.filter[key].includes(elem));
      }
    )
  ) {
    return prevQuery;
  }
  prevQuery = query;
  return query;
}

function toEntity(item: RawOntologyTerm) {
  const [id, name] = Object.entries(item)[0];

  return { id, name: name || id || "" };
}

function toCitationEntity(citation: string) {
  return {
    id: citation,
    name: formatCitation(citation),
  };
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

interface MarkerGeneRequestBody {
  celltype: string;
  n_markers: number;
  organism: string;
  test: "ttest" | "binomtest";
  tissue: string;
}

interface HardcodedMarkerGeneRequest extends MarkerGeneRequestBody {
  n_markers: 25;
}

export function generateMarkerGeneBody(
  cellTypeID: string,
  tissueID: string,
  organismID: string,
  test: "ttest" | "binomtest"
): HardcodedMarkerGeneRequest {
  return {
    celltype: cellTypeID,
    n_markers: 25,
    organism: organismID,
    test: test,
    tissue: tissueID,
  };
}

export interface FetchMarkerGeneParams {
  cellTypeID: string;
  tissueID: string;
  organismID: string;
  test?: "ttest" | "binomtest";
}

export interface FetchGeneInfoParams {
  geneID: string;
  geneSymbol: string;
  signal?: AbortSignal;
}

export async function fetchMarkerGenes({
  cellTypeID,
  organismID,
  tissueID,
  test = "ttest",
}: FetchMarkerGeneParams): Promise<MarkerGeneResponse> {
  const url = API_URL + API.WMG_MARKER_GENES;
  const body = generateMarkerGeneBody(cellTypeID, tissueID, organismID, test);
  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
    body: JSON.stringify(body),
    method: "POST",
  });

  const json: MarkerGeneResponse = await response.json();

  if (!response.ok) {
    throw json;
  }

  return json;
}

export async function fetchGeneInfo({
  geneID,
  geneSymbol,
  signal,
}: FetchGeneInfoParams): Promise<GeneInfo> {
  const url =
    API_URL + API.WMG_GENE_INFO + `?geneID=${geneID}&gene=${geneSymbol}`;
  const response = await fetch(url, {
    ...DEFAULT_FETCH_OPTIONS, // Required for CORS
    signal,
  });

  const json: GeneInfo = await response.json();

  if (!response.ok) {
    throw json;
  }
  return json;
}

export const USE_MARKER_GENES = {
  ENTITIES: [ENTITIES.WMG_MARKER_GENES],
  id: "wmg-marker-genes",
};

export const USE_GENE_INFO = {
  ENTITIES: [ENTITIES.WMG_GENE_INFO],
  id: "wmg-gene-info",
};

export interface MarkerGenesByCellType {
  [cellType: string]: MarkerGeneResponse["marker_genes"];
}

export interface MarkerGene {
  gene_ontology_term_id: string;
  effect_size: number;
  p_value: number;
}

export interface MarkerGeneResponse<T = MarkerGene[]> {
  marker_genes: T;
  snapshot_id: string;
}

export function useMarkerGenes({
  cellTypeID,
  organismID,
  tissueID,
  test,
}: FetchMarkerGeneParams): UseQueryResult<MarkerGeneResponse<MarkerGene>> {
  const { data } = usePrimaryFilterDimensions(2);
  const genesByID = useMemo((): { [name: string]: OntologyTerm } => {
    let result: { [name: string]: OntologyTerm } = {};

    if (!data) return result;

    const { genes } = data;

    result = generateTermsByKey(genes, "id");

    return result;
  }, [data]);

  function filterMarkerGenes(markerGenes: MarkerGene[]): MarkerGene[] {
    return markerGenes.filter(
      (markerGene) => markerGene.effect_size >= FMG_GENE_STRENGTH_THRESHOLD
    );
  }

  return useQuery(
    /**
     * (thuang): Add all arguments to `fetchMarkerGenes()` as dependencies,
     * so React Query can cache responses correctly without running into
     * issues like #4161
     */
    [USE_MARKER_GENES, cellTypeID, organismID, test, tissueID, 2],
    async () => {
      const output = await fetchMarkerGenes({
        cellTypeID,
        organismID,
        test,
        tissueID,
      });
      const markerGenesIndexedByGeneName = Object.fromEntries(
        filterMarkerGenes(output.marker_genes).reduce(
          (newEntries, { gene_ontology_term_id, ...data }) => {
            try {
              newEntries.push([genesByID[gene_ontology_term_id].name, data]);
            } catch (e) {
              console.log("could not find gene with id", gene_ontology_term_id);
            }
            return newEntries;
          },
          [] as [string, Partial<MarkerGene>][]
        )
      );
      return { ...output, marker_genes: markerGenesIndexedByGeneName };
    },
    {
      staleTime: Infinity,
    }
  );
}

export function useGeneInfo(
  geneSymbol: string,
  version: 1 | 2 = 2
): UseQueryResult<GeneInfo> {
  const { data } = usePrimaryFilterDimensions(version);
  const genesByName = useMemo((): { [name: string]: OntologyTerm } => {
    let result: { [name: string]: OntologyTerm } = {};

    if (!data) return result;

    const { genes } = data;

    result = generateTermsByKey(genes, "name");

    return result;
  }, [data]);

  const geneID = genesByName[geneSymbol]?.id;
  return useQuery(
    // Sometimes using the UUID doesn't work if it's out-of-date and so the gene symbol is a fallback.
    [USE_GENE_INFO, geneSymbol, geneID],
    ({ signal }) => fetchGeneInfo({ geneSymbol, geneID, signal }),
    {
      staleTime: Infinity,
    }
  );
}

export const CELL_TYPE_VIEW_ID_DIVIDER = "$";

function getCellTypeViewId(
  ontologyTermId: CellTypeId,
  optionId: CompareOptionId
): ViewId {
  return `${ontologyTermId}${CELL_TYPE_VIEW_ID_DIVIDER}${optionId}`;
}

export function getOntologyTermIdFromCellTypeViewId(
  viewId: ViewId
): CompareOptionId {
  return viewId.split(CELL_TYPE_VIEW_ID_DIVIDER)[0];
}

export function getOptionIdFromCellTypeViewId(viewId: ViewId): CompareOptionId {
  return viewId.split(CELL_TYPE_VIEW_ID_DIVIDER)[1];
}
