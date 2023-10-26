import { useCallback, useContext, useEffect, useMemo, useState } from "react";
import { Organism as IOrganism } from "src/views/WheresMyGene/common/types";
import { useAvailableOrganisms } from "src/common/queries/wheresMyGene";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import { selectOrganism } from "src/views/WheresMyGene/common/store/actions";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { selectHasCustomFiltersOrGenesSelected } from "src/views/WheresMyGene/common/store/selectors";

export const useConnect = () => {
  const dispatch = useContext(DispatchContext);
  const state = useContext(StateContext);
  const { selectedOrganismId } = state;

  const { data: organisms } = useAvailableOrganisms(2);

  const [isDialogOpen, setIsDialogOpen] = useState(false);

  const [pendingOrganism, setPendingOrganism] = useState<IOrganism | null>(
    null
  );

  const hasCustomFiltersOrGenesSelected =
    selectHasCustomFiltersOrGenesSelected(state);

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

    setPendingOrganism(organism);
  }

  const handleDialogConfirm = useCallback(() => {
    track(EVENTS.WMG_SELECT_ORGANISM, { payload: pendingOrganism?.name });

    dispatch?.(selectOrganism(pendingOrganism?.id || null));

    setIsDialogOpen(false);
    setPendingOrganism(null);
  }, [dispatch, pendingOrganism]);

  function handleDialogCancel() {
    setPendingOrganism(null);
    setIsDialogOpen(false);
  }

  useEffect(() => {
    if (!pendingOrganism) return;

    if (hasCustomFiltersOrGenesSelected) {
      setIsDialogOpen(true);
    } else {
      handleDialogConfirm();
    }
  }, [pendingOrganism, hasCustomFiltersOrGenesSelected, handleDialogConfirm]);

  const organism = organismsById[selectedOrganismId || ""];

  /**
   * FIXME(thuang): We need to pass `value` to Dropdown this way instead of using
   * `Dropdown.props.value`, because when we cancel the confirm dialog, the
   * internal Dropdown selected value is still staying as the new value instead of
   * reverting back to the old value, even though the old value is still being
   * passed to Dropdown.
   * To fix this, SDS Dropdown needs to detect that it's in a controlled state,
   * if so, it should use the value passed to it instead of the internal state.
   */
  const DropdownMenuProps = useMemo(
    () => ({
      value: organism,
    }),
    [organism]
  );

  return {
    handleDialogCancel,
    handleDialogConfirm,
    handleOnChange,
    isDialogOpen,
    organism,
    organisms,
    setIsDialogOpen,
    DropdownMenuProps,
  };
};
