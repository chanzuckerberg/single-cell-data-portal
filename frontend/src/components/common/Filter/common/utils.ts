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
 * Returns formatted number with corresponding magnitude.
 * @param num - Number to format.
 * @returns String containing number with corresponding magnitude.
 */
export function formatNumberToMagnitude(num: number): string {
  // TODO(cc) demo code only - must be productionalized.
  const magnitude = num < 100000 ? 1000 : 1000000;
  const precision = num < 1000000 ? 1 : 2;
  const symbol = magnitude === 1000 ? "k" : "M";
  return `${roundToPrecision(precision, num) / magnitude}${symbol}`;
}

/**
 * Round the given number to the given significant digits.
 * @param precision - Number of significant digits.
 * @param num - Number to round to significant digits.
 * @returns Rounded number to the given significant digits.
 */
function roundToPrecision(precision: number, num: number): number {
  return parseFloat(num.toPrecision(precision));
}
