/**
 * Helper functions for the controls reducer
 */

import {
  CategoricalAnnotationColumnSchema,
  Category,
  Schema,
} from "../../common/types/schema";
import { DataframeColumn } from "../dataframe";

import fromEntries from "../fromEntries";
import { isCategoricalAnnotation } from "./annotationsHelpers";

/*
Selection state for categoricals are tracked in an Object that
has two main components for each category:
1. mapping of option value to an index
2. array of bool selection state by index
Remember that option values can be ANY js type, except undefined/null.

  {
    _category_name_1: {
      // map of option value to index
      categoryValueIndices: Map([
        catval1: index,
        ...
      ])

      // index->selection true/false state
      categoryValueSelected: [ true/false, true/false, ... ]

      // number of options
      numCategoryValues: number,
    }
  }
*/

export function isSelectableCategoryName(
  schema: Schema,
  name: string,
): boolean {
  const { index } = schema.annotations.obs;
  const colSchema = schema.annotations.obsByName[name];
  return Boolean(
    name &&
      name !== index &&
      (isCategoricalAnnotation(schema, name) || colSchema.writable),
  );
}

/**
 * Returns all obs annotation names that are categorical AND have a
 * "reasonably" small number of categories AND are not the index column.
 * If the initial name list not provided, use everything in the schema.
 */
export function selectableCategoryNames(
  schema: Schema,
  names?: string[],
): string[] {
  if (!schema) return [];

  if (!names) names = schema.annotations.obs.columns.map((c) => c.name);

  return names.filter((name) => isSelectableCategoryName(schema, name));
}

export interface CategorySummary {
  // Array of natively typed category values (all of them)
  allCategoryValues: Category[];
  // Array of natively typed category values (top N only)
  categoryValues: Category[];
  // Category value (native type) -> category index (top N only)
  categoryValueIndices: Map<Category, number>;
  // Number of values in the category (top N)
  numCategoryValues: number;
  // Array of cardinality of each category, (top N)
  categoryValueCounts: number[];
  isUserAnno: boolean;
}

/**
 * Summarize the annotation data currently in dataframe column.  Must return
 * categoryValues in sorted order, and must include all category values even
 * if they are not actively used in the current annoMatrix view.
 */
export function createCategorySummaryFromDfCol(
  dfCol: DataframeColumn,
  colSchema: CategoricalAnnotationColumnSchema,
): CategorySummary {
  const { writable: isUserAnno } = colSchema;
  const summary = dfCol.summarizeCategorical();
  const { categories: allCategoryValues } = colSchema;
  const categoryValues = allCategoryValues;
  const categoryValueCounts = allCategoryValues.map(
    (cat: Category) => summary.categoryCounts.get(cat) ?? 0,
  );
  const categoryValueIndices = new Map(categoryValues.map((v, i) => [v, i]));

  const numCategoryValues = categoryValueIndices.size;

  return {
    allCategoryValues,
    categoryValues,
    categoryValueIndices,
    numCategoryValues,
    categoryValueCounts,
    isUserAnno,
  };
}

/**
 * key is label, value is selected or not
 * See client/src/reducers/categoricalSelection.ts for state description
 */
type SelectionState = Map<string, boolean>;

export interface CategoricalSelection {
  // Categorical annotation name
  [key: string]: SelectionState;
}

export function createCategoricalSelection(
  names: string[],
): CategoricalSelection {
  return fromEntries(names.map((name) => [name, new Map()]));
}
