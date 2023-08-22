import { useMemo } from "react";
import { EnrichedGenesQueryResponse } from "src/common/queries/cellGuide";
import { FMG_GENE_STRENGTH_THRESHOLD } from "src/views/WheresMyGene/common/constants";
import { HOMO_SAPIENS } from "../constants";

interface ComputationalMarkerGeneTableData {
  symbol: string;
  name: string;
  marker_score: string;
  me: string;
  pc: string;
}

function _getSortedOrganisms(genes: EnrichedGenesQueryResponse): string[] {
  const organisms = new Set<string>();
  for (const markerGene of genes) {
    if (markerGene.marker_score < FMG_GENE_STRENGTH_THRESHOLD) continue;
    organisms.add(markerGene.organism);
  }
  return Array.from(organisms).sort((a, b) => {
    if (a === HOMO_SAPIENS) return -1;
    if (b === HOMO_SAPIENS) return 1;
    return a.localeCompare(b);
  });
}

function _getSortedOrgans(
  genes: EnrichedGenesQueryResponse,
  selectedOrganismFilter: string
): string[] {
  const organs = new Set<string>();
  for (const markerGene of genes) {
    if (markerGene.organism !== selectedOrganismFilter) continue;
    if (markerGene.marker_score < FMG_GENE_STRENGTH_THRESHOLD) continue;
    organs.add(markerGene.tissue);
  }

  return Array.from(organs).sort((a, b) => {
    return a.localeCompare(b);
  });
}

export function useComputationalMarkerGenesTableRowsAndFilters({
  genes,
  selectedOrganism,
  selectedOrgan,
}: {
  genes: EnrichedGenesQueryResponse;
  selectedOrganism: string;
  selectedOrgan: string;
}): {
  computationalMarkerGeneTableData: ComputationalMarkerGeneTableData[];
} {
  return useMemo(() => {
    if (!genes)
      return {
        computationalMarkerGeneTableData: [],
      };

    // get sorted organisms
    const sortedOrganisms = _getSortedOrganisms(genes);
    const selectedOrganismFilter = !sortedOrganisms.includes(selectedOrganism)
      ? ""
      : selectedOrganism;

    // get sorted organs
    const sortedOrgans = _getSortedOrgans(genes, selectedOrganismFilter);
    const selectedOrganFilter = !sortedOrgans.includes(selectedOrgan)
      ? ""
      : selectedOrgan;

    const rows: ComputationalMarkerGeneTableData[] = [];
    for (const markerGene of genes) {
      const { pc, me, name, symbol, organism, marker_score, tissue } =
        markerGene;
      if (organism !== selectedOrganismFilter) continue;
      if (tissue !== selectedOrganFilter) continue;
      if (marker_score < FMG_GENE_STRENGTH_THRESHOLD) continue;
      rows.push({
        symbol,
        name,
        marker_score: marker_score.toFixed(2),
        me: me.toFixed(2),
        pc: (pc * 100).toFixed(1),
      });
    }

    return {
      computationalMarkerGeneTableData: rows,
    };
  }, [genes, selectedOrganism, selectedOrgan]);
}
