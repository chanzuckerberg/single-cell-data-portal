import { useCallback, useContext } from "react";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGeneV2/common/store";
import {
  resetGenesToDelete,
  toggleGeneToDelete,
} from "src/views/WheresMyGeneV2/common/store/actions";

export function useDeleteGenes(): {
  genesToDelete: string[];
  handleGeneClick: (gene: string) => void;
  reset: () => void;
} {
  const dispatch = useContext(DispatchContext);
  const { genesToDelete } = useContext(StateContext);

  const handleGeneClick = useCallback(
    (gene: string) => {
      if (!dispatch) return;

      dispatch(toggleGeneToDelete(gene));
    },
    [dispatch]
  );

  const reset = useCallback(() => {
    if (!dispatch) return;

    dispatch(resetGenesToDelete());
  }, [dispatch]);

  return {
    genesToDelete,
    handleGeneClick,
    reset,
  };
}
