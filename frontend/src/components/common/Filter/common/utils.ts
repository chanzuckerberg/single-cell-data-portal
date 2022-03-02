// App dependencies
import { Ontology } from "src/common/entities";
import {
  Categories,
  OntologyCategoryKey,
} from "src/components/common/Filter/common/entities";

/**
 * Return type of accessor function.
 */
type OntologyCellAccessorFn = (categories: Categories) => string[];

/**
 * Number of digits to format decmial points to.
 */
const FIXED_TO = 1;

/**
 * Value to divide numbers by when scaling by million.
 */
const SCALE_MILLION = 1000000;

/**
 * Value to divide numbers by when scaling by thousand.
 */
const SCALE_THOUSAND = 1000;

/**
 * Scale symbol for million.
 */
export const SYMBOL_MILLION = "M";

/**
 * Scale symbol for million.
 */
export const SYMBOL_THOUSAND = "k";

/**
 * Create function to be used by column.accessor in react-table column definition, for columns containing ontology
 * metadata (ontology label and key) values.
 * @param categoryKey - Object key of value to display in cell.
 * @returns Function that returns value with the given key, to display in a cell.
 */
export function ontologyCellAccessorFn(
  categoryKey: OntologyCategoryKey
): OntologyCellAccessorFn {
  return (categories: Categories) =>
    categories[categoryKey].map((o: Ontology) => o.label);
}

/**
 * Format number to the correct scale and possibly round.
 * @param num - Number to format.
 * @returns String representation of given number, possibly rounded, with corresponding scale symbol.
 */
export function formatNumberToScale(num: number): string {
  // Return numbers less than 1000 rounded.
  if (num < 1000) {
    return `${Math.round(num)}`;
  }

  // Handle numbers smaller than 10,000: format to thousands. For example, 1400 becomes 1.4k.
  if (num < 10000) {
    const formatted = scaleAndFix(num, SCALE_THOUSAND, FIXED_TO);
    return `${formatted}${SYMBOL_THOUSAND}`;
  }

  // Handle numbers smaller than 1,000,000: format to thousands and round. For example, 14,500 becomes 15k.
  if (num < 1000000) {
    const rounded = Math.round(scaleAndFix(num, SCALE_THOUSAND, FIXED_TO));
    return `${rounded}${SYMBOL_THOUSAND}`;
  }

  // Handle numbers larger than 1,000,000: format to millions. For example, 1,450,000 becomes 1.5M.
  const formatted = scaleAndFix(num, SCALE_MILLION, FIXED_TO);
  return `${formatted}${SYMBOL_MILLION}`;
}

/**
 * Round up the given number to the given fixed digits, removing any insignificant decimals.
 * @param num - Number to round to significant digits.
 * @param scale - Number to divide given number by.
 * @param fixTo - Number of digits to fix number to.
 * @returns Rounded number to the given fixed digits.
 */
function scaleAndFix(num: number, scale: number, fixTo: number): number {
  return parseFloat((num / scale + Number.EPSILON).toFixed(fixTo));
}
