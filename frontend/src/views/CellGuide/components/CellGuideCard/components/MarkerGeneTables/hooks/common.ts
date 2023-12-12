import { useMemo } from "react";
import {
  useComputationalMarkers,
  useAllOrgansLookupTables,
} from "src/common/queries/cellGuide";
import { ALL_TISSUES, HOMO_SAPIENS, NO_ORGAN_ID } from "../constants";

/* NOTE: Here "organ" refers to a coarse tissue or system (ex: "eye", "digestive system")
   whereas, "tissue" refers coarse tissue/system AND granular tissue. That is, "tissue"
   refers to a superset that contains "organs"
   (ex: "retina", "lacrimal gland" are granular tissues that belong to the "eye" organ)
*/

/* NOTE: Currently uses ONLY computational marker genes to derive the
   filter list of organs and organisms. In the future, we might use an
   independent filter list of organs and organisms.
*/
export function useOrganAndOrganismFilterListForCellType(cellTypeId: string): {
  organsMap: Map<string, string>;
  organismsList: string[];
  isSuccess: boolean;
} {
  const {
    data: computationalMarkers,
    isSuccess: isComputationalMarkersSuccess,
  } = useComputationalMarkers(cellTypeId);

  const { data: organLabelToIdMap, isSuccess } = useAllOrgansLookupTables();

  // eslint-disable-next-line sonarjs/cognitive-complexity
  return useMemo(() => {
    const organisms = new Set<string>([HOMO_SAPIENS]);
    const filteredOrgansMap = new Map<string, string>([
      [ALL_TISSUES, NO_ORGAN_ID],
    ]);
    let organId: string | undefined;

    // 1. construct a label to id map of the organs intersected with
    //    the computational marker genes.
    // 2. construct a list of unique organisms in the enriched genes.
    if (computationalMarkers) {
      for (const markerGene of computationalMarkers) {
        const organLabel = markerGene.groupby_dims.tissue_ontology_term_label;
        if (organLabel && organLabelToIdMap) {
          organId = organLabelToIdMap.get(organLabel);
          if (organId) {
            filteredOrgansMap.set(organLabel, organId);
          }
        }
        organisms.add(markerGene.groupby_dims.organism_ontology_term_label);
      }
    }

    // sort the organs map derived from the enriched genes
    const sortedFilteredOrganMap = new Map(
      [...filteredOrgansMap.entries()].sort((a, b) => {
        if (a[0] === ALL_TISSUES) return -1;
        if (b[0] === ALL_TISSUES) return 1;
        return a[0].localeCompare(b[0]);
      })
    );

    // sort the organisms list derived from the enriched genes
    const sortedOrganismList = Array.from(organisms).sort((a, b) => {
      if (a === HOMO_SAPIENS) return -1;
      if (b === HOMO_SAPIENS) return 1;
      return a.localeCompare(b);
    });

    return {
      organsMap: sortedFilteredOrganMap,
      organismsList: sortedOrganismList,
      /**
       * (thuang): Expose `isSuccess`, so `CellGuide/components/CellGuideCard/connect.ts`
       * can use it to determine if the data is ready and determine if the user should
       * be redirected to the tissue agnostic cell type page.
       */
      isSuccess: isComputationalMarkersSuccess && isSuccess,
    };
  }, [
    computationalMarkers,
    organLabelToIdMap,
    isSuccess,
    isComputationalMarkersSuccess,
  ]);
}
