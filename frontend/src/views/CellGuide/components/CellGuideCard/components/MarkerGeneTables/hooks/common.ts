import { useMemo } from "react";
import { useEnrichedGenes } from "src/common/queries/cellGuide";
import { ALL_TISSUES, HOMO_SAPIENS } from "../constants";

export function useMarkerGenesTableTissueAndOrganismFilterListForCelltype(
  cellTypeId: string
): {
  uniqueTissues: string[];
  uniqueOrganisms: string[];
} {
  const { data: enrichedGenes } = useEnrichedGenes(cellTypeId);

  return useMemo(() => {
    const tissues = new Set<string>([ALL_TISSUES]);
    const organisms = new Set<string>([HOMO_SAPIENS]);

    if (enrichedGenes) {
      for (const markerGene of enrichedGenes) {
        tissues.add(markerGene.tissue);
        organisms.add(markerGene.organism);
      }
    }

    const tissueList = Array.from(tissues).sort((a, b) => {
      if (a === ALL_TISSUES) return -1;
      if (b === ALL_TISSUES) return 1;
      return a.localeCompare(b);
    });

    const organismList = Array.from(organisms).sort((a, b) => {
      if (a === HOMO_SAPIENS) return -1;
      if (b === HOMO_SAPIENS) return 1;
      return a.localeCompare(b);
    });

    return {
      uniqueTissues: tissueList,
      uniqueOrganisms: organismList,
    };
  }, [enrichedGenes]);
}
