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
 * Based off https://stackoverflow.com/a/9462382.
 * @param num
 * @returns string (a displayable number with corresponding magnitude).
 */
export function formatNumberToMagnitude(num: number): string {
  const si = [
    { symbol: "", value: 1 },
    { symbol: "k", value: 1e3 },
    { symbol: "M", value: 1e6 },
    { symbol: "G", value: 1e9 },
    { symbol: "T", value: 1e12 },
    { symbol: "PB", value: 1e15 },
    { symbol: "E", value: 1e18 },
  ];
  const rx = /\.0+$|(\.[0-9]*[1-9])0+$/;
  let i;
  for (i = si.length - 1; i > 0; i--) {
    if (num >= si[i].value) {
      break;
    }
  }
  const symbol = si[i].symbol;
  return (num / si[i].value).toFixed(1).replace(rx, "$1") + symbol;
}
