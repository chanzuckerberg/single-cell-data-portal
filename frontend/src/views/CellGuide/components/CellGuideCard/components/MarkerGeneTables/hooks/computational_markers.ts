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

enum selectionCriteriaErrorCode {
  OrganismMismatch,
  TissueMismatch,
  LowMarkerScore,
  Success,
}

// Apply selection criteria in the following order and short circut when a criteria fails:
// 1. Organism
// 2. Tissue
// 3. Marker score
function _applyOrderedSelectionCriteria({
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
}): { pass: boolean; errorCode: selectionCriteriaErrorCode } {
  const { groupby_dims, marker_score } = markerGene;

  // When groupby_dims.tissue_ontology_term_label is undefined, default to ALL_TISSUES
  const {
    organism_ontology_term_label,
    tissue_ontology_term_label = ALL_TISSUES,
  } = groupby_dims;

  // 1. Filter by organism
  if (organism_ontology_term_label !== selectedOrganismLabel)
    return {
      pass: false,
      errorCode: selectionCriteriaErrorCode.OrganismMismatch,
    };

  // 2. Filter by tissue
  //
  // There are marker genes tissues labeled as "All Tissues" AND there
  // there are marker genes tissues that are undefined which are re-labeled
  // as ALL_TISSUES. Select them only when selectedOrganLabel is "All Tissues"
  let pass = false;

  if (selectedOrganLabel === ALL_TISSUES) {
    pass = tissue_ontology_term_label === ALL_TISSUES;
  } else {
    const tissueId = allTissuesLabelToIdMap.get(tissue_ontology_term_label);
    pass = Boolean(
      tissueId &&
        isTissueIdDescendantOfAncestorTissueId(tissueId, selectedOrganId)
    );
  }

  if (!pass)
    return {
      pass: false,
      errorCode: selectionCriteriaErrorCode.TissueMismatch,
    };

  // 3. Finally, Filter by marker score
  if (marker_score < FMG_GENE_STRENGTH_THRESHOLD)
    return {
      pass: false,
      errorCode: selectionCriteriaErrorCode.LowMarkerScore,
    };

  return {
    pass: true,
    errorCode: selectionCriteriaErrorCode.Success,
  };
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
  allFilteredByLowMarkerScore: boolean;
} {
  return useMemo(() => {
    if (!genes)
      return {
        computationalMarkerGeneTableData: [],
        allFilteredByLowMarkerScore: false,
      };

    const rows: ComputationalMarkerGeneTableData[] = [];
    let allFilteredByLowMarkerScore = false;

    for (const markerGene of genes) {
      const { pass, errorCode } = _applyOrderedSelectionCriteria({
        markerGene: markerGene,
        allTissuesLabelToIdMap: allTissuesLabelToIdMap,
        selectedOrganismLabel: selectedOrganismLabel,
        selectedOrganLabel: selectedOrganLabel,
        selectedOrganId: selectedOrganId,
      });

      if (!pass) {
        // Because of the ordering of the selection criteria, we can identify whether
        // the set of genes that pass the organism and tissue filter but entirely fail
        // the marker score filter.
        //
        // ignore genes that that fail the organism or tissue filter
        if (
          errorCode === selectionCriteriaErrorCode.OrganismMismatch ||
          errorCode === selectionCriteriaErrorCode.TissueMismatch
        ) {
          continue;
        }

        // of those genes that pass the organism and tissue filter,
        // check if they all fail the marker score filter.
        // This check is done for future proofing - it is possible we
        // may introduce other filter checks in the future.
        allFilteredByLowMarkerScore =
          errorCode === selectionCriteriaErrorCode.LowMarkerScore;
        continue;
      }

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
      allFilteredByLowMarkerScore: allFilteredByLowMarkerScore,
    };
  }, [
    genes,
    allTissuesLabelToIdMap,
    selectedOrganismLabel,
    selectedOrganLabel,
    selectedOrganId,
  ]);
}
