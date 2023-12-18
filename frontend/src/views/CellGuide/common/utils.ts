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

import { ROUTES } from "src/common/constants/routes";
import { NO_ORGAN_ID } from "src/views/CellGuide/components/CellGuideCard/components/MarkerGeneTables/constants";

export function getCellTypeLink({
  /**
   * (thuang): `queryTissueId` can be undefined when a user is on `/cellguide/:cellTypeId` route
   * instead of `/cellguide/tissues/:tissueId/cell-types/:cellTypeId` route.
   *
   * NOTE: `tissue` and `organ` in variable names are interchangeable here.
   */
  tissueId = NO_ORGAN_ID,
  cellTypeId,
}: {
  tissueId: string;
  cellTypeId: string;
}) {
  const urlCellTypeId = cellTypeId.replace(":", "_") ?? "";
  const urlTissueId = tissueId.replace(":", "_") || NO_ORGAN_ID;

  if (tissueId === NO_ORGAN_ID) {
    return ROUTES.CELL_GUIDE_CELL_TYPE.replace(":cellTypeId", urlCellTypeId);
  } else {
    return ROUTES.CELL_GUIDE_TISSUE_SPECIFIC_CELL_TYPE.replace(
      ":tissueId",
      urlTissueId
    ).replace(":cellTypeId", urlCellTypeId);
  }
}
