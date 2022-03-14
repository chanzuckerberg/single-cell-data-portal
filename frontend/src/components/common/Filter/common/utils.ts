// App dependencies
import { Ontology } from "src/common/entities";
import {
  Categories,
  OntologyCategoryKey,
  OntologyNode,
  SPECIES_KEY,
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
 * //TODO(cc)
 * Find the IDs of all ancestors of the given ontology.
 * @param ontologyNode - Node to of ontology to find ancestors of.
 * @param ancestors - The set of ancestors of the given ontology ID.
 * @returns Array of ontology IDs of nodes that are ancestors of the given ontology ID.
 */
export function findOntologyAncestorIds(
  rootNodes: OntologyNode[],
  ontologyNode: OntologyNode,
  ancestors: string[]
): OntologyNode | undefined {
  for (let i = 0; i < rootNodes.length; i++) {
    const rootNode = rootNodes[i];
    if (rootNode === ontologyNode) {
      return ontologyNode;
    }
    if (!rootNode.children) {
      continue;
    }
    let descendantNode;
    for (let j = 0; j < rootNode.children.length; j++) {
      const childNode = rootNode.children[j];
      descendantNode = findOntologyAncestorIds(
        [childNode],
        ontologyNode,
        ancestors
      );
      if (descendantNode) {
        break;
      }
    }
    if (descendantNode) {
      ancestors.push(rootNode.ontology_term_id);
      return descendantNode;
    }
  }
}

/**
 * Find the IDs of all descendants of the given ontology.
 * @param ontologyNode - Node to of ontology to find descendants of.
 * @param descendants - The set of descendants of the given ontology ID.
 * @returns Array of ontology IDs of nodes that are descendants of the given ontology ID.
 */
export function findOntologyDescendantIds(
  ontologyNode: OntologyNode,
  descendants: string[] = []
): string[] {
  if (!ontologyNode.children) {
    return descendants;
  }

  ontologyNode.children.forEach((childNode) => {
    descendants.push(childNode.ontology_term_id);
    findOntologyDescendantIds(childNode, descendants);
  });

  return descendants;
}

/**
 * Determine the ontology key of the given ontology ID. For example, "HsapDv:0000003" returns "HsapDv".
 * @param ontologyId - ID to determine ontology key from.
 * @returns String containing ontology key.
 */
export function getOntologySpeciesKey(ontologyId: string): SPECIES_KEY {
  return SPECIES_KEY[ontologyId.split(":")[0] as keyof typeof SPECIES_KEY];
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
