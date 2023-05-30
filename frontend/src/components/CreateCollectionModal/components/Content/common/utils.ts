/**
 * Utilities for collection form-related functionality.
 */

import { COLLECTION_LINK_TYPE } from "src/common/entities";
import { DefaultDropdownMenuOption } from "@czi-sds/components";
import { Value as DropdownValue } from "src/components/common/Form/Dropdown";
import { CONSORTIA } from "src/components/CreateCollectionModal/components/Content/common/constants";

/**
 * Returns consortia, reformatted in a suitable shape for dropdown menu select options.
 * @returns consortia menu select options.
 */
export function buildConsortiaOptions(
  consortia: CONSORTIA[]
): DefaultDropdownMenuOption[] {
  return consortia.map((consortium) => {
    return {
      name: consortium,
    };
  });
}

/**
 * Determine if given link type is DOI.
 * @param linkType - Type of link to check if DOI.
 * @returns True if link type is DOI.
 */
export function isLinkTypeDOI(linkType: COLLECTION_LINK_TYPE): boolean {
  return linkType === COLLECTION_LINK_TYPE.DOI;
}

/**
 * Sorts selected consortia options alphabetically.
 * @param consortia - Selected consortia options.
 * @returns selected consorted options sorted alphabetically.
 */
export function sortConsortia(
  consortia: DropdownValue
): DefaultDropdownMenuOption[] {
  if (Array.isArray(consortia)) {
    return consortia.sort(({ name: c0 }, { name: c1 }) => {
      if (c0 > c1) {
        return 1;
      }
      if (c0 < c1) {
        return -1;
      }
      return 0;
    });
  }
  return [];
}
