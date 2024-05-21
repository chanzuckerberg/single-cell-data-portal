import { useCallback, useContext, useEffect, useMemo } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { useAvailableOrganisms } from "src/common/queries/differentialExpression";
import {
  DispatchContext,
  StateContext,
} from "src/views/DifferentialExpression/common/store";
import {
  clearQueryGroup1Filters,
  clearQueryGroup2Filters,
  clearSubmittedQueryGroups,
  selectOrganism,
} from "src/views/DifferentialExpression/common/store/actions";
import { Organism as IOrganism } from "src/views/DifferentialExpression/common/types";

export const useConnect = () => {
  const dispatch = useContext(DispatchContext);
  const { organismId } = useContext(StateContext);
  const { data: organisms } = useAvailableOrganisms();

  const handleOnChange = useCallback(
    (organism: IOrganism | null): void => {
      if (!dispatch || !organism || organismId === organism.id) return;

      track(EVENTS.DE_SELECT_ORGANISM, { payload: organism?.name });
      dispatch(clearQueryGroup1Filters());
      dispatch(clearQueryGroup2Filters());
      dispatch(clearSubmittedQueryGroups());
      dispatch(selectOrganism(organism?.id || null));
    },
    [dispatch, organismId]
  );

  useEffect(() => {
    if (!organisms || !dispatch || organismId) return;

    const organism = organisms.find(
      (organism: IOrganism) => organism.name === "Homo sapiens"
    );

    if (!organism) return;

    dispatch(selectOrganism(organism.id));
  }, [organisms, dispatch, organismId]);

  const organismsById = useMemo(() => {
    const result: { [id: string]: IOrganism } = {};

    if (!organisms) return result;

    for (const organism of organisms) {
      result[organism.id] = organism;
    }

    return result;
  }, [organisms]);

  const organism = organismsById[organismId || ""];

  return {
    handleOnChange,
    organism,
    organisms,
  };
};
