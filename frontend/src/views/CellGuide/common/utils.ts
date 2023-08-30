import { TISSUE_DESCENDANTS } from "src/components/common/Filter/common/constants";

export function isTissueIdDescendantOfAncestorTissueId(
  tissueId: string,
  ancestorTissueId: string
): boolean {
  // a tissue is a descendant of itself
  if (tissueId === ancestorTissueId) return true;

  // check all descendants of ancestorTissueId excluding ancestorTissueId itself
  const descendantTermIds = TISSUE_DESCENDANTS[ancestorTissueId];
  if (!descendantTermIds) return false;

  return descendantTermIds.includes(tissueId);
}

export function filterDescendantsOfAncestorTissueId(
  tissueIdList: string[],
  ancestorTissueId: string
): string[] {
  const descendantList = [ancestorTissueId];

  const moreDescendantTermIds = TISSUE_DESCENDANTS[ancestorTissueId];
  if (moreDescendantTermIds) {
    descendantList.push(...moreDescendantTermIds);
  }
  const descendantSet = new Set(descendantList);

  return tissueIdList.filter((value) => descendantSet.has(value));
}
