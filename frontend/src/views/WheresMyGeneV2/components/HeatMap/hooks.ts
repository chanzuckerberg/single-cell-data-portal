import { Dispatch, useEffect, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { autoExpandTissues } from "src/views/WheresMyGene/common/store/actions";
import { PayloadAction } from "src/views/WheresMyGene/common/store/reducer";

interface Props {
  selectedGenes: string[];
  selectedCellTypes: string[];
  displayedCellTypes: Set<string>;
}

export function useTrackHeatMapLoaded({
  selectedGenes,
  selectedCellTypes,
  displayedCellTypes,
}: Props): void {
  /**
   * (thuang): We only want to send `WMG_HEATMAP_LOADED` event when at least
   * one gene are selected, and if the flag `isHeatMapEventFired` is `false`.
   * We reset `isHeatMapEventFired` to `false` when the user removes all the
   * selected tissues and/or all the genes.
   *
   * For example:
   * 1. Visitor visits WMG
   * 2. Visitor selects at least one gene
   * 3. Heat map is loaded --> triggers one heat map event. `isHeatMapEventFired` is set to `true`
   * 4. Visitor adds another tissue or gene, nothing happens
   * 5. Visitor removes all tissues, heat map disappears. `isHeatMapEventFired` is set to `false`
   * 6. Visitor then adds another tissue and the gene is still selected
   * 7. Another heat map appears --> triggers one heat map event. `isHeatMapEventFired` is set to `true`
   * 8. Visitor removes all genes, heat map shows only the tissue names. `isHeatMapEventFired` is set to `false`
   * 9. Visitor then adds another gene --> triggers one heat map event. `isHeatMapEventFired` is set to `true`
   * Total: 3 heat map events
   */
  useEffect(() => {
    const hasSelectedGenes = selectedGenes.length > 0;

    if (!hasSelectedGenes) return;

    track(EVENTS.WMG_HEATMAP_LOADED, {
      num_cell_types_selected: selectedCellTypes.length,
      num_rows_include: displayedCellTypes.size,
    });
  }, [selectedGenes, selectedCellTypes, displayedCellTypes]);
}

/**
 * (thuang): WMG Reducer handles most of the `expandedTissueIds` logic, but the
 * following cases are not handled by the reducer, since only the Heatmap component
 * can detect the following states
 */
export function useHandleExpandedTissueIds({
  filteredCellTypeIds,
  filteredTissueIds,
  displayedTissueIds,
  dispatch,
}: {
  filteredCellTypeIds: string[];
  filteredTissueIds: string[];
  displayedTissueIds: string[];
  dispatch: Dispatch<PayloadAction<unknown>> | null;
}) {
  const [prevFilteredCellTypeIds, setPrevFilteredCellTypeIds] =
    useState(filteredCellTypeIds);

  const [prevFilteredTissueIds, setPrevFilteredTissueIds] =
    useState(filteredTissueIds);

  // DEBUG
  // DEBUG
  // DEBUG
  console.log("-👀---------displayedTissueIds", displayedTissueIds);
  console.log("-👀---------prevFilteredCellTypeIds", prevFilteredCellTypeIds);
  console.log("-👀---------filteredCellTypeIds", filteredCellTypeIds);
  console.log("-👀---------prevFilteredTissueIds", prevFilteredTissueIds);
  console.log("-👀---------filteredTissueIds", filteredTissueIds);

  useEffect(() => {
    /**
     * (thuang): When we exit the tissue filter mode, but still have cell type filter,
     * we want to expand all the tissues.
     */
    if (
      prevFilteredTissueIds.length &&
      !filteredTissueIds.length &&
      filteredCellTypeIds.length
    ) {
      expandTissues();
    }

    /**
     * (thuang): When we exit the cell type filter mode, but still have tissue filter,
     * we want to expand all the tissues.
     */
    if (
      prevFilteredCellTypeIds.length &&
      !filteredCellTypeIds.length &&
      filteredTissueIds.length
    ) {
      expandTissues();
    }

    if (prevFilteredCellTypeIds !== filteredCellTypeIds) {
      setPrevFilteredCellTypeIds(filteredCellTypeIds);
    }

    if (prevFilteredTissueIds !== filteredTissueIds) {
      setPrevFilteredTissueIds(filteredTissueIds);
    }

    function expandTissues() {
      dispatch?.(autoExpandTissues(displayedTissueIds));
    }
  }, [
    filteredCellTypeIds,
    prevFilteredCellTypeIds,
    filteredTissueIds,
    prevFilteredTissueIds,
    displayedTissueIds,
    dispatch,
  ]);
}
