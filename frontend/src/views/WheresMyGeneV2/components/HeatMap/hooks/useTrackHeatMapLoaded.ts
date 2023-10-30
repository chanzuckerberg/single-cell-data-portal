import { useEffect, useRef } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";

interface Props {
  selectedTissues: string[];
  selectedGenes: string[];
}

export function useTrackHeatMapLoaded({
  selectedTissues,
  selectedGenes,
}: Props): void {
  const isHeatMapEventFired = useRef(false);

  /**
   * (thuang): We only want to send `WMG_HEATMAP_LOADED` event when at least one
   * tissue and one gene are selected, and if the flag `isHeatMapEventFired` is `false`.
   * We reset `isHeatMapEventFired` to `false` when the user removes all the
   * selected tissues and/or all the genes.
   *
   * For example:
   * 1. Visitor visits WMG
   * 2. Visitor selects at least one tissue and at least one gene
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
    const hasSelectedTissues = selectedTissues.length > 0;
    const hasSelectedGenes = selectedGenes.length > 0;

    if (!hasSelectedTissues || !hasSelectedGenes) {
      isHeatMapEventFired.current = false;
      return;
    }

    if (
      !isHeatMapEventFired.current &&
      hasSelectedTissues &&
      hasSelectedGenes
    ) {
      track(EVENTS.WMG_HEATMAP_LOADED);
      isHeatMapEventFired.current = true;
    }
  }, [selectedTissues, selectedGenes]);
}
