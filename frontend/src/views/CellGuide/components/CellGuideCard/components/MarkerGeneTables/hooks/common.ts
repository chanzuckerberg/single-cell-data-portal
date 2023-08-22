import { useMemo, useState } from "react";
import {
  CanonicalMarkersQueryResponse,
  EnrichedGenesQueryResponse,
  useCanonicalMarkers,
  useEnrichedGenes,
} from "src/common/queries/cellGuide";
import { ALL_TISSUES, HOMO_SAPIENS } from "../constants";

export function useMarkerGenesTableTissueAndOrganismFilterListForCelltype(
  cellTypeId: string
): {
  uniqueTissues: string[];
  uniqueOrganisms: string[];
} {
  const [computationalMarkerGenes, setComputationalMarkerGenes] =
    useState<EnrichedGenesQueryResponse>([]);
  const [canonicalMarkerGenes, setCanonicalMarkerGenes] =
    useState<CanonicalMarkersQueryResponse>([]);

  const { data: enrichedGenes } = useEnrichedGenes(cellTypeId);
  const { data: canonicalMarkers } = useCanonicalMarkers(cellTypeId);

  return useMemo(() => {
    const tissues = new Set<string>([ALL_TISSUES]);
    const organisms = new Set<string>([HOMO_SAPIENS]);

    if (canonicalMarkers) {
      setCanonicalMarkerGenes(canonicalMarkers);
    }

    if (enrichedGenes) {
      setComputationalMarkerGenes(enrichedGenes);
    }

    for (const markerGene of canonicalMarkerGenes) {
      tissues.add(markerGene.tissue);
    }

    for (const markerGene of computationalMarkerGenes) {
      tissues.add(markerGene.tissue);
      organisms.add(markerGene.organism);
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
  }, [
    canonicalMarkerGenes,
    canonicalMarkers,
    computationalMarkerGenes,
    enrichedGenes,
  ]);
}
