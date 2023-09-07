import { CATEGORY_FILTER_ID } from "src/components/common/Filter/common/entities";
import { VIEW_MODE } from "src/common/hooks/useViewMode";

/**
 * Collection cell count object key.
 */
export const COLLECTION_CELL_COUNT = "cell_count";

/**
 * Collection curator name object key.
 */
export const COLLECTION_CURATOR_NAME = "curator_name";

/**
 * Collection ID object key.
 */
export const COLLECTION_ID = "id";

/**
 * Collection name object key.
 */
export const COLLECTION_NAME = "name";

/**
 * Collection recency object key.
 */
export const COLLECTION_RECENCY = "recency";

/**
 * Collection revised by object key.
 */
export const COLLECTION_REVISED_BY = "revisedBy";

/**
 * Collection summary citation object key.
 */
export const COLLECTION_SUMMARY_CITATION = "summaryCitation";

/**
 * Collection status object key.
 */
export const COLLECTION_STATUS = "status";

/**
 * Key identifying recency sort by column.
 */
export const COLUMN_ID_RECENCY = "recency";

/**
 * Collections column deny list.
 */
export const COLLECTIONS_COLUMN_DENY_LIST: (string | CATEGORY_FILTER_ID)[] = [
  COLLECTION_ID,
  COLUMN_ID_RECENCY,
  COLLECTION_REVISED_BY,
  CATEGORY_FILTER_ID.CELL_TYPE_CALCULATED,
  CATEGORY_FILTER_ID.DEVELOPMENT_STAGE,
  CATEGORY_FILTER_ID.PUBLICATION,
  CATEGORY_FILTER_ID.PUBLICATION_DATE_VALUES,
  CATEGORY_FILTER_ID.SELF_REPORTED_ETHNICITY,
  CATEGORY_FILTER_ID.SEX,
  CATEGORY_FILTER_ID.SUSPENSION_TYPE,
  CATEGORY_FILTER_ID.TISSUE_CALCULATED,
];

/**
 * Category filter id deny list for collections mode.
 */
export const COLLECTIONS_MODE_CATEGORY_FILTER_DENY_LIST: CATEGORY_FILTER_ID[] =
  [
    CATEGORY_FILTER_ID.CELL_COUNT,
    CATEGORY_FILTER_ID.CURATOR_NAME,
    CATEGORY_FILTER_ID.GENE_COUNT,
    CATEGORY_FILTER_ID.STATUS,
  ];

/**
 * Column deny list for collections mode.
 */
export const COLLECTIONS_MODE_COLUMN_DENY_LIST: (
  | string
  | CATEGORY_FILTER_ID
)[] = [
  ...COLLECTIONS_COLUMN_DENY_LIST,
  CATEGORY_FILTER_ID.ASSAY,
  CATEGORY_FILTER_ID.CELL_COUNT,
  CATEGORY_FILTER_ID.CURATOR_NAME,
  CATEGORY_FILTER_ID.STATUS,
];

/**
 * Category filter id deny list for curator mode.
 */
export const CURATOR_MODE_CATEGORY_FILTER_DENY_LIST: CATEGORY_FILTER_ID[] = [
  CATEGORY_FILTER_ID.CELL_COUNT,
  CATEGORY_FILTER_ID.GENE_COUNT,
];

/**
 * Category filter id partition list for curator mode.
 */
export const CURATOR_MODE_CATEGORY_FILTER_PARTITION_LIST: CATEGORY_FILTER_ID[] =
  [CATEGORY_FILTER_ID.CURATOR_NAME, CATEGORY_FILTER_ID.STATUS];

/**
 * Column deny list for curator mode.
 */
export const CURATOR_MODE_COLUMN_DENY_LIST: (string | CATEGORY_FILTER_ID)[] = [
  ...COLLECTIONS_COLUMN_DENY_LIST,
  COLLECTION_SUMMARY_CITATION,
];

/**
 * Category filter id deny list.
 */
export const CATEGORY_FILTER_DENY_LIST: Record<
  VIEW_MODE,
  CATEGORY_FILTER_ID[]
> = {
  [VIEW_MODE.DEFAULT]: COLLECTIONS_MODE_CATEGORY_FILTER_DENY_LIST,
  [VIEW_MODE.CURATOR]: CURATOR_MODE_CATEGORY_FILTER_DENY_LIST,
};

/**
 * Category filter id partition list.
 */
export const CATEGORY_FILTER_PARTITION_LIST: Record<
  VIEW_MODE,
  CATEGORY_FILTER_ID[]
> = {
  [VIEW_MODE.DEFAULT]: [],
  [VIEW_MODE.CURATOR]: CURATOR_MODE_CATEGORY_FILTER_PARTITION_LIST,
};

/**
 * Column deny list.
 */
export const COLUMN_DENY_LIST: Record<
  VIEW_MODE,
  (string | CATEGORY_FILTER_ID)[]
> = {
  [VIEW_MODE.DEFAULT]: COLLECTIONS_MODE_COLUMN_DENY_LIST,
  [VIEW_MODE.CURATOR]: CURATOR_MODE_COLUMN_DENY_LIST,
};
