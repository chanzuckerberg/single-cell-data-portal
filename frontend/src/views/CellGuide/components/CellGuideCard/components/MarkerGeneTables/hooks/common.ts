import { useMemo } from "react";
import {
  useAllTissuesLabelToIdMap,
  useEnrichedGenes,
} from "src/common/queries/cellGuide";
import { ALL_TISSUES, HOMO_SAPIENS, NO_TISSUE_ID } from "../constants";

export function useMarkerGenesTableTissueAndOrganismFilterListForCelltype(
  cellTypeId: string
): {
  tissuesMap: Map<string, string>;
  organismsList: string[];
} {
  const { data: enrichedGenes } = useEnrichedGenes(cellTypeId);

  // get tissue-label: tissue-id map for all tissues
  const allTissuesMap = useAllTissuesLabelToIdMap();

  return useMemo(() => {
    const organisms = new Set<string>([HOMO_SAPIENS]);
    const filteredTissuesMap = new Map<string, string>([
      [ALL_TISSUES, NO_TISSUE_ID],
    ]);
    let tissueId: string | undefined;

    // 1. construct a tissue-label: tissue-id map of the tissues in the enriched genes
    // 2. construct a list of unique organisms in the enriched genes
    if (enrichedGenes) {
      for (const markerGene of enrichedGenes) {
        tissueId = allTissuesMap.get(markerGene.tissue);
        if (tissueId) {
          filteredTissuesMap.set(markerGene.tissue, tissueId);
        }
        organisms.add(markerGene.organism);
      }
    }

    // sort the tissues map derived from the enriched genes
    const sortedFilteredTissueMap = new Map(
      [...filteredTissuesMap.entries()].sort((a, b) => {
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
      tissuesMap: sortedFilteredTissueMap,
      organismsList: sortedOrganismList,
    };
  }, [enrichedGenes, allTissuesMap]);
}
