import { useMemo } from "react";
import {
  ComputationalMarkersQueryResponse,
  ComputationalMarkersQueryResponseEntry,
} from "src/common/queries/cellGuide";
import { FMG_GENE_STRENGTH_THRESHOLD } from "src/views/WheresMyGene/common/constants";
import { isTissueIdDescendantOfAncestorTissueId } from "src/views/CellGuide/common/utils";
import { ALL_TISSUES } from "../constants";

interface ComputationalMarkerGeneTableData {
  symbol: string;
  name: string;
  marker_score: string;
  me: string;
  pc: string;
}

function _passSelectionCriteria({
  markerGene,
  allTissuesLabelToIdMap,
  selectedOrganismLabel,
  selectedOrganLabel,
  selectedOrganId,
}: {
  markerGene: ComputationalMarkersQueryResponseEntry;
  selectedOrganismLabel: string;
  selectedOrganLabel: string;
  selectedOrganId: string;
  allTissuesLabelToIdMap: Map<string, string>;
}): boolean {
  const { groupby_dims, marker_score } = markerGene;

  // When groupby_dims.tissue_ontology_term_label is undefined, default to ALL_TISSUES
  const {
    organism_ontology_term_label,
    tissue_ontology_term_label = ALL_TISSUES,
  } = groupby_dims;

  if (marker_score < FMG_GENE_STRENGTH_THRESHOLD) return false;
  if (organism_ontology_term_label !== selectedOrganismLabel) return false;

  // filter by tissue
  //
  // There are marker genes tissues labeled as "All Tissues" AND there
  // there are marker genes tissues that are undefined which are re-labeled
  // as ALL_TISSUES. Select them only when selectedOrganLabel is "All Tissues"
  if (selectedOrganLabel === ALL_TISSUES) {
    return tissue_ontology_term_label === ALL_TISSUES;
  }

  const tissue_id = allTissuesLabelToIdMap.get(tissue_ontology_term_label);

  if (!tissue_id) return false;

  return isTissueIdDescendantOfAncestorTissueId(tissue_id, selectedOrganId);
}

export function useComputationalMarkerGenesTableRowsAndFilters({
  genes,
  allTissuesLabelToIdMap,
  selectedOrganismLabel,
  selectedOrganLabel,
  selectedOrganId,
}: {
  genes: ComputationalMarkersQueryResponse;
  allTissuesLabelToIdMap: Map<string, string>;
  selectedOrganismLabel: string;
  selectedOrganLabel: string;
  selectedOrganId: string;
}): {
  computationalMarkerGeneTableData: ComputationalMarkerGeneTableData[];
} {
  return useMemo(() => {
    if (!genes)
      return {
        computationalMarkerGeneTableData: [],
      };

    const rows: ComputationalMarkerGeneTableData[] = [];

    for (const markerGene of genes) {
      if (
        !_passSelectionCriteria({
          markerGene: markerGene,
          allTissuesLabelToIdMap: allTissuesLabelToIdMap,
          selectedOrganismLabel: selectedOrganismLabel,
          selectedOrganLabel: selectedOrganLabel,
          selectedOrganId: selectedOrganId,
        })
      )
        continue;

      const { pc, me, name, symbol, marker_score } = markerGene;

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
  }, [
    genes,
    allTissuesLabelToIdMap,
    selectedOrganismLabel,
    selectedOrganLabel,
    selectedOrganId,
  ]);
}
