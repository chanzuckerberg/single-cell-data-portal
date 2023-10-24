import { useContext, useEffect, useMemo } from "react";
import { Organism as IOrganism } from "src/views/WheresMyGeneV2/common/types";
import { useAvailableOrganisms } from "src/common/queries/wheresMyGene";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGeneV2/common/store";
import { selectOrganism } from "src/views/WheresMyGeneV2/common/store/actions";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";

export const useConnect = () => {
  const dispatch = useContext(DispatchContext);
  const { selectedOrganismId } = useContext(StateContext);

  const { data: organisms } = useAvailableOrganisms(2);

  // (thuang): Default to "Homo sapiens" on first load
  useEffect(() => {
    if (!organisms || !dispatch || selectedOrganismId) return;

    const organism = organisms.find(
      (organism: IOrganism) => organism.name === "Homo sapiens"
    );

    if (!organism) return;

    dispatch(selectOrganism(organism.id));
  }, [organisms, dispatch, selectedOrganismId]);

  const organismsById = useMemo(() => {
    const result: { [id: string]: IOrganism } = {};

    if (!organisms) return result;

    for (const organism of organisms) {
      result[organism.id] = organism;
    }

    return result;
  }, [organisms]);

  function handleOnChange(organism: IOrganism | null): void {
    if (!dispatch || !organism || selectedOrganismId === organism.id) return;

    track(EVENTS.WMG_SELECT_ORGANISM, { payload: organism?.name });

    dispatch(selectOrganism(organism?.id || null));
  }

  const organism = organismsById[selectedOrganismId || ""];
  return { organisms, organism, handleOnChange };
};
