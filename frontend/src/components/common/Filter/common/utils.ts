// App dependencies
import { Ontology } from "src/common/entities";
import {
  Categories,
  CategoriesKeyOfTypeOntologyArray,
  ONTOLOGY_VIEW_KEY,
  OntologyNode,
  OntologyTermSet,
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
 * Find the node with the given ontology ID.
 * @param rootNodes - Top-level nodes in ontology tree to search for ontology ID in.
 * @param ontologyId - The ontology ID of the node to find.
 */
export function findOntologyNodeById(
  rootNodes: OntologyNode[],
  ontologyId: string
): OntologyNode | undefined {
  for (let i = 0; i < rootNodes.length; i++) {
    const rootNode = rootNodes[i];
    if (rootNode.ontology_term_id === ontologyId) {
      return rootNode;
    }
    const node = findOntologyNodeById(rootNode.children ?? [], ontologyId);
    if (node) {
      return node;
    }
  }
}

/**
 * Find the parent node of the given node.
 * @param rootNodes - Node tree structure to find parent in.
 * @param ontologyNode - Node to find parent of.
 * @returns Ontology node that is the parent of the given ontology node.
 */
export function findOntologyParentNode(
  rootNodes: OntologyNode[],
  ontologyNode: OntologyNode
): OntologyNode | undefined {
  for (let i = 0; i < rootNodes.length; i++) {
    // Check if this node is the parent.
    const rootNode = rootNodes[i];
    const isParent = rootNode.children?.find(
      (childNode) => childNode === ontologyNode
    );
    if (isParent) {
      return rootNode;
    }

    // Otherwise, check each child to see if they are the parent.
    const parentNode = findOntologyParentNode(
      rootNode.children ?? [],
      ontologyNode
    );
    if (parentNode) {
      return parentNode;
    }
  }
}

/**
 * Determine the ontology key of the given ontology ID. For example, "HsapDv:0000003" returns "HsapDv".
 * @param ontologyId - ID to determine ontology key from.
 * @returns String containing ontology key.
 */
export function getOntologySpeciesKey(ontologyId: string): ONTOLOGY_VIEW_KEY {
  return ONTOLOGY_VIEW_KEY[
    ontologyId.split(":")[0] as keyof typeof ONTOLOGY_VIEW_KEY
  ];
}

/**
 * List all ontology IDs in the given ontology node.
 * @param ontologyIds - Set to add ontology IDs of node to.
 * @param node - Ontology node to list ontology IDs of.
 */
function listOntologyNodeIds(ontologyIds: Set<string>, node: OntologyNode) {
  // Add this node to the set.
  ontologyIds.add(node.ontology_term_id);

  // Find ontology IDs for each child of this node.
  node.children?.forEach((childNode) =>
    listOntologyNodeIds(ontologyIds, childNode)
  );
}

/**
 * List all ontology IDs in the given ontology tree.
 * @param ontologyTermSet - Ontology view model to list ontology IDs of.
 * @returns Set of all ontology IDs present in the given ontology tree.
 */
export function listOntologyTreeIds(
  ontologyTermSet: OntologyTermSet
): Set<string> {
  return Object.values(ontologyTermSet).reduce((accum, node) => {
    node.forEach((node) => listOntologyNodeIds(accum, node));
    return accum;
  }, new Set<string>());
}

/**
 * Create function to be used by column.accessor in react-table column definition, for columns containing ontology
 * metadata (ontology label) values.
 * @param categoryKey - Object key of value to display in cell.
 * @returns Function that returns label of ontology value with the given key, to display in a cell.
 */
export function ontologyLabelCellAccessorFn<
  K extends CategoriesKeyOfTypeOntologyArray
>(categoryKey: K): OntologyCellAccessorFn {
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
