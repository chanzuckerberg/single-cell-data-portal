/**
 * Utilities for collection form-related functionality.
 */

import { COLLECTION_LINK_TYPE } from "src/common/entities";

/**
 * Determine if given link type is DOI.
 * @param linkType - Type of link to check if DOI.
 * @returns True if link type is DOI.
 */
export function isLinkTypeDOI(linkType: COLLECTION_LINK_TYPE): boolean {
  return linkType === COLLECTION_LINK_TYPE.DOI;
}
