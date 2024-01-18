import { useCallback, useState } from "react";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { FEATURES } from "src/common/featureFlags/features";

export type OnReorderFn = (mode: REORDER_MODE) => void;

/**
 * Reorder mode.
 */
export enum REORDER_MODE {
  ACTIVE = "ACTIVE",
  INACTIVE = "INACTIVE",
}

interface UseReorderMode {
  isReorderUX: boolean;
  mode: REORDER_MODE;
  onReorder: OnReorderFn;
}

/**
 * Determines the reorder datasets mode for the collection.
 * The reorder mode can be either "inactive" or "active" and is used to enable or disable the datasets reorder feature
 * in the collection view.
 * @returns reorder mode.
 */
export function useReorderMode(): UseReorderMode {
  const isReorderUX = useFeatureFlag(FEATURES.REORDER); // Reorder datasets UX feature flag (reordering is currently only available with the feature flag).
  const [mode, setMode] = useState<REORDER_MODE>(REORDER_MODE.INACTIVE);

  // Updates reorder mode.
  const onReorder = useCallback((mode: REORDER_MODE) => {
    setMode(mode);
  }, []);

  return { isReorderUX, mode, onReorder };
}
