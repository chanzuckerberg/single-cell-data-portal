// App dependencies
import { Ontology } from "src/common/entities";
import {
  Categories,
  DevelopmentStageNode,
  DEVELOPMENT_STAGES,
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
const SCALE_MILLION = 1_000_000;

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
 * Development stage leaf node ontology IDs keyed by each ancestor node ontology ID, used to facilitate ancestor-leaf
 * lookups.
 */
export const DEVELOPMENT_STAGE_LEAF_ONTOLOGY_IDS_BY_ANCESTOR = Object.values(
  DEVELOPMENT_STAGES
).reduce((accum, developmentStages) => {
  developmentStages.forEach((developmentStage) => {
    setLeafOntologyIdsForDevelopmentStage(accum, developmentStage);
  });
  return accum;
}, new Map<string, Set<string>>());

/**
 * Set of all development stage ontology IDs, used by filter functionality to determine the full set of possible
 * ontology IDs.
 * TODO(cc) share recursion with setLeafOntologyIdsForDevelopmentStage but without passing in accum?
 * TODO(cc) test
 */
export const DEVELOPMENT_STAGE_ONTOLOGY_IDS = [
  ...DEVELOPMENT_STAGE_LEAF_ONTOLOGY_IDS_BY_ANCESTOR.keys(),
].reduce((accum, ontologyId) => {
  accum.add(ontologyId);
  DEVELOPMENT_STAGE_LEAF_ONTOLOGY_IDS_BY_ANCESTOR.get(ontologyId)?.forEach(
    (childOntologyId) => {
      accum.add(childOntologyId);
    }
  );
  return accum;
}, new Set<string>());

/**
 * Find all leaf node (ontology IDs) for the given development stage node.
 * @param developmentStage - Development stage node to find leaf node ontology IDs of.
 * @param leafOntologyIds - Set of leaf node ontology IDs found for the given development stage.
 */
function findLeafNodeOntologyIds(
  developmentStage: DevelopmentStageNode,
  leafOntologyIds: Set<string>
) {
  // If this development stage node has no children, we have found a leaf! Add ontology ID to set.
  if (!developmentStage.children) {
    leafOntologyIds.add(developmentStage.ontology_term_id);
    return;
  }
  developmentStage.children.forEach((childDevelopmentStage) =>
    findLeafNodeOntologyIds(childDevelopmentStage, leafOntologyIds)
  );
}

/**
 * Format number to the correct scale and possibly round.
 * @param num - Number to format.
 * @returns String representation of given number, possibly rounded, with corresponding scale symbol.
 */
export function formatNumberToScale(num: number): string {
  // Return numbers less than 1000 rounded.
  if (num < SCALE_THOUSAND) {
    return `${Math.round(num)}`;
  }

  // Handle numbers smaller than 10,000: format to thousands. For example, 1400 becomes 1.4k.
  if (num < 10 * SCALE_THOUSAND) {
    const formatted = scaleAndFix(num, SCALE_THOUSAND, FIXED_TO);
    return `${formatted}${SYMBOL_THOUSAND}`;
  }

  // Handle numbers smaller than 1,000,000: format to thousands and round. For example, 14,500 becomes 15k.
  if (num < SCALE_MILLION) {
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

/**
 * Find the leaf nodes for non-leaf development stage nodes in the tree of the given node and add to map.
 * @param leafOntologyIdsByAncestor - Map of leaf ontology IDs keyed by ancestor development stage ontology ID.
 * @param developmentStage - Development stage node to find leaf node ontology IDs of.
 */
function setLeafOntologyIdsForDevelopmentStage(
  leafOntologyIdsByAncestor: Map<string, Set<string>>,
  developmentStage: DevelopmentStageNode
) {
  // Node itself a leaf.
  if (!developmentStage.children) {
    return;
  }

  // Find the leaf nodes for this development stage node.
  const leafNodes = new Set<string>();
  findLeafNodeOntologyIds(developmentStage, leafNodes);
  leafOntologyIdsByAncestor.set(developmentStage.ontology_term_id, leafNodes);

  // Find the leaf nodes for children of this development stage node, if any.
  developmentStage.children.forEach((childDevelopmentStage) =>
    setLeafOntologyIdsForDevelopmentStage(
      leafOntologyIdsByAncestor,
      childDevelopmentStage
    )
  );
}
