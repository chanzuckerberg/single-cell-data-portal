import { useMemo } from "react";
import { useAllOrgansLookupTables } from "src/common/queries/cellGuide";
import { ALL_TISSUES, NO_ORGAN_ID } from "../constants";

/* NOTE: Here "organ" refers to a coarse tissue or system (ex: "eye", "digestive system")
   whereas, "tissue" refers coarse tissue/system AND granular tissue. That is, "tissue"
   refers to a superset that contains "organs"
   (ex: "retina", "lacrimal gland" are granular tissues that belong to the "eye" organ)
*/

/* NOTE: Currently uses ONLY computational marker genes to derive the
   filter list of organs and organisms. In the future, we might use an
   independent filter list of organs and organisms.
*/
export function useOrganAndOrganismFilterListForCellType(): {
  organsMap: Map<string, string>;
  organismsList: string[];
  isSuccess: boolean;
} {
  const { data: organLabelToIdMap, isSuccess } = useAllOrgansLookupTables();
  // eslint-disable-next-line sonarjs/cognitive-complexity
  return useMemo(() => {
    organLabelToIdMap.set(ALL_TISSUES, NO_ORGAN_ID);
    const sortedFilteredOrganMap = new Map(
      [...organLabelToIdMap.entries()].sort((a, b) => {
        if (a[0] === ALL_TISSUES) return -1;
        if (b[0] === ALL_TISSUES) return 1;
        return a[0].localeCompare(b[0]);
      })
    );

    return {
      organsMap: sortedFilteredOrganMap,
      // hardcode this to the only two organisms present in Census for now
      // TODO: get this from the backend
      organismsList: ["Homo sapiens", "Mus musculus"],
      /**
       * (thuang): Expose `isSuccess`, so `CellGuide/components/CellGuideCard/connect.ts`
       * can use it to determine if the data is ready and determine if the user should
       * be redirected to the tissue agnostic cell type page.
       */
      isSuccess,
    };
  }, [organLabelToIdMap, isSuccess]);
}
