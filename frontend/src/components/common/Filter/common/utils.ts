// App dependencies
import { Ontology } from "src/common/entities";
import {
  Categories,
  OntologyCategoryKey,
} from "src/components/common/Filter/common/entities";

/* Return type of accessor function. */
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
