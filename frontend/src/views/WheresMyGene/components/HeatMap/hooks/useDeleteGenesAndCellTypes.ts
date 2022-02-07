import { useCallback, useContext } from "react";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import {
  resetGenesToDeleteAndCellTypeIdsToDelete,
  toggleCellTypeIdToDelete,
  toggleGeneToDelete,
} from "src/views/WheresMyGene/common/store/actions";

export function useDeleteGenesAndCellTypes(): {
  cellTypeIdsToDelete: string[];
  genesToDelete: string[];
  handleCellTypeClick: (cellTypeId: string) => void;
  handleGeneClick: (gene: string) => void;
  reset: () => void;
} {
  const dispatch = useContext(DispatchContext);
  const { genesToDelete, cellTypeIdsToDelete } = useContext(StateContext);

  const handleGeneClick = useCallback(
    (gene: string) => {
      if (!dispatch) return;

      dispatch(toggleGeneToDelete(gene));
    },
    [dispatch]
  );

  const handleCellTypeClick = useCallback(
    (cellTypeId: string) => {
      if (!dispatch) return;

      dispatch(toggleCellTypeIdToDelete(cellTypeId));
    },
    [dispatch]
  );

  const reset = useCallback(() => {
    if (!dispatch) return;

    dispatch(resetGenesToDeleteAndCellTypeIdsToDelete());
  }, [dispatch]);

  return {
    cellTypeIdsToDelete,
    genesToDelete,
    handleCellTypeClick,
    handleGeneClick,
    reset,
  };
}
