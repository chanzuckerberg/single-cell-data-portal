import { useMemo } from "react";
import {
  useComputationalMarkers,
  useTissueMetadata,
} from "src/common/queries/cellGuide";
import { ALL_TISSUES, HOMO_SAPIENS, NO_ORGAN_ID } from "../constants";

/* Currently uses ONLY computational marker genes (ie enriched genes)
   to derive the filter list of tissues and organisms
*/
export function useOrganAndOrganismFilterListForCelltype(cellTypeId: string): {
  organsMap: Map<string, string>;
  organismsList: string[];
} {
  const { data: computationalMarkers } = useComputationalMarkers(cellTypeId);

  // get label: ontology-term-id map for all tissues
  // only a subset of the organs in this map will be used
  // to construct the tissue filter list
  const { data: allOrgansMap } = useTissueMetadata();
  // eslint-disable-next-line sonarjs/cognitive-complexity
  return useMemo(() => {
    const organisms = new Set<string>([HOMO_SAPIENS]);
    const filteredOrgansMap = new Map<string, string>([
      [ALL_TISSUES, NO_ORGAN_ID],
    ]);

    // 1. construct a label: id map of the tissues in the enriched genes
    // 2. construct a list of unique organisms in the enriched genes
    if (computationalMarkers) {
      for (const markerGene of computationalMarkers) {
        organisms.add(markerGene.groupby_dims.organism_ontology_term_label);
      }
    }

    if (allOrgansMap) {
      for (const tissueId in allOrgansMap) {
        filteredOrgansMap.set(allOrgansMap[tissueId].name, tissueId);
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
    };
  }, [computationalMarkers, allOrgansMap]);
}
