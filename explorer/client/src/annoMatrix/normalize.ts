import { _getColumnSchema, _isIndex } from "./schema";
import catLabelSort from "../util/catLabelSort";
import {
  unassignedCategoryLabel,
  overflowCategoryLabel,
  globalConfig,
} from "../globals";
import { Dataframe, LabelType, DataframeColumn } from "../util/dataframe";
import {
  AnnotationColumnSchema,
  Field,
  Schema,
  Category,
  CategoricalAnnotationColumnSchema,
} from "../common/types/schema";
import { isDictEncodedTypedArray } from "../common/types/arraytypes";

export function normalizeResponse(
  field: Field,
  schema: Schema,
  response: Dataframe,
): Dataframe {
  /**
   * There are a number of assumptions in the front-end about data typing and data
   * characteristics. This routine will normalize a server response dataframe
   * to match front-end expectations and UI conventions. This includes cast/transform
   * of the data and schema updates.
   *
   * This consolidates all assumptions into one location, for ease of update.
   *
   * Currently, this includes normalization for obs/var columns only:
   *
   * - Dataframe columns in var/obs that are declared type: boolean may be sent by
   *   the server in a variety of formats (eg uint8, etc).  Cast to JS Array[boolean]
   *
   * - "Categorical" columns may not have all categories represented in the server-provided
   *   schema (for valid reasons, eg, floating point rounding differences).  For all
   *   types we treat as categorical in the UI (string, boolean, categorical), update
   *   the schema to contain all categories as a convenience.
   *
   * - "Categorical" columns (ie, string, boolean, categorical) may contain an excess
   *   of category values (aka labels).  Consolidate any excess into an "all other"
   *   category.
   */

  // currently no data or schema normalization necessary for X or emb
  if (field !== Field.obs && field !== Field.var) return response;

  const colLabels = response.colIndex.labels();
  for (const colLabel of colLabels) {
    const colSchema = _getColumnSchema(
      schema,
      field,
      colLabel,
    ) as AnnotationColumnSchema;

    const isIndex = _isIndex(schema, field, colLabel);
    const { type, writable } = colSchema;

    // Boolean data -- cast entire array to Array[bool]
    if (type === "boolean") {
      response = castColumnToBoolean(response, colLabel);
    }

    // Types that are categorical in UI (string, boolean, categorical) OR are writable
    // are introspected to ensure the schema `categories` field and data values match,
    // and that we do not have an excess of category values (for non-writable columns)
    const isEnumType =
      type === "boolean" ||
      type === "string" ||
      type === "categorical" ||
      writable;
    // If "feature_id" column exists, should not normalize
    // (gene info cards require ensembl ID of any gene, and normalizing overwrites some of those IDs)
    if (!isIndex && colLabel !== "feature_id" && isEnumType) {
      response = normalizeCategorical(response, colLabel, colSchema);
    }
  }
  return response;
}

function castColumnToBoolean(df: Dataframe, label: LabelType): Dataframe {
  const colData = df.col(label).asArray();
  const newColData = new Array(colData.length);
  for (let i = 0; i < colData.length; i += 1) newColData[i] = !!colData[i];
  df = df.replaceColData(label, newColData);
  return df;
}

export function normalizeWritableCategoricalSchema(
  colSchema: AnnotationColumnSchema,
  col: DataframeColumn,
): CategoricalAnnotationColumnSchema {
  /*
  Ensure all enum writable / categorical schema have a categories array, that
  the categories array contains all unique values in the data array, AND that
  the array is UI sorted.
  */
  const categorySet = new Set<Category>(
    col.summarizeCategorical().categories.concat(
      // TODO #35: Use type guards instead of casting
      (colSchema as CategoricalAnnotationColumnSchema).categories ?? [],
    ),
  );
  if (!categorySet.has(unassignedCategoryLabel)) {
    categorySet.add(unassignedCategoryLabel);
  }
  // TODO #35: Use type guards instead of casting

  (colSchema as CategoricalAnnotationColumnSchema).categories = catLabelSort(
    true,
    Array.from(categorySet),
  );
  return colSchema as CategoricalAnnotationColumnSchema;
}

export function normalizeCategorical(
  df: Dataframe,
  colLabel: LabelType,
  colSchema: AnnotationColumnSchema,
): Dataframe {
  /*
  If writable, ensure schema matches data and we have an unassigned label

  If not writable, ensure schema matches data and that we consolidate labels in excess
  of "top N" into an overflow labels.
  */
  const { writable } = colSchema;
  const col = df.col(colLabel);
  if (writable) {
    // writable (aka user) annotations
    normalizeWritableCategoricalSchema(colSchema, col);
    return df;
  }

  // else read-only, categorical columns
  const TopN = globalConfig.maxCategoricalOptionsToDisplay;

  // consolidate all categories from data and schema into a single list
  const colDataSummary = col.summarizeCategorical();
  const allCategories = new Set<Category>(colDataSummary.categories);

  // if no overflow, just UI sort schema categories and return
  if (allCategories.size <= TopN) {
    (colSchema as CategoricalAnnotationColumnSchema).categories = catLabelSort(
      writable,
      [...allCategories.keys()],
    );
    return df;
  }

  // Otherwise, pick top N categories by count and rewrite data

  // choose unique overflow category label
  let overflowCatName = `${colLabel}${overflowCategoryLabel}`;
  while (allCategories.has(overflowCatName)) {
    overflowCatName += "_";
  }

  // pick top N category labels and add overflow label
  const topNCategories = new Set(
    [...colDataSummary.categoryCounts.keys()].slice(0, TopN),
  );
  topNCategories.add(overflowCatName);

  // rewrite data - consolidate all excess labels into overflow label
  const colData = col.asArray();
  const newColData = new Array(colData.length);
  for (let i = 0; i < newColData.length; i += 1) {
    if (!topNCategories.has(colData[i])) {
      newColData[i] = overflowCatName;
    } else if (isDictEncodedTypedArray(colData)) {
      newColData[i] = colData.vat(i);
    } else {
      newColData[i] = colData[i];
    }
  }
  // replace data in dataframe
  df = df.replaceColData(colLabel, newColData);

  // Update schema with categories, in UI sort order. Ensure overflow label is at end
  // of list for display purposes.
  const revisedCategories = df.col(colLabel).summarizeCategorical().categories;
  revisedCategories.push(
    revisedCategories.splice(revisedCategories.indexOf(overflowCatName), 1)[0],
  );
  // TODO #35: Use type guards instead of casting
  (colSchema as CategoricalAnnotationColumnSchema).categories = catLabelSort(
    writable,
    revisedCategories,
  );

  return df;
}
