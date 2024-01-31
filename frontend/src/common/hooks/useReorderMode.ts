import { useCallback, useEffect, useState } from "react";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { FEATURES } from "src/common/featureFlags/features";

export type OnCancelReorderFn = () => void;

export type OnReorderFn = (
  datasetID: string,
  targetDatasetID: string,
  position: ORDER_POSITION
) => void;

export type OnSaveReorderFn = () => void;

export type OnStartReorderFn = () => void;

export enum ORDER_POSITION {
  BEFORE = -1,
  AFTER = 1,
}

export interface ReorderAction {
  onCancelReorder: OnCancelReorderFn;
  onReorder: OnReorderFn;
  onSaveReorder: OnSaveReorderFn;
  onStartReorder: OnStartReorderFn;
}

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
  orderedIDs: string[];
  reorderAction: ReorderAction;
}

/**
 * Reorder functionality for collection datasets.
 * The reorder mode can be either "inactive" or "active" and is used to enable or disable the datasets reorder feature
 * in the collection view.
 * @param initialOrderedIDs - Initial dataset IDs, ordered.
 * @returns reorder mode.
 */
export function useReorderMode(initialOrderedIDs?: string[]): UseReorderMode {
  const isReorderUX = useFeatureFlag(FEATURES.REORDER); // Reorder datasets UX feature flag (reordering is currently only available with the feature flag).
  const [mode, setMode] = useState<REORDER_MODE>(REORDER_MODE.INACTIVE);
  const [orderedIDs, setOrderedIDs] = useState<string[]>([]);

  // Cancels reorder mode.
  const onCancelReorder = useCallback(() => {
    setMode(REORDER_MODE.INACTIVE);
    setOrderedIDs(initialOrderedIDs || []);
  }, [initialOrderedIDs]);

  // Updates order.
  const onReorder = useCallback(
    (
      datasetID: string,
      targetDatasetID: string,
      orderPosition: ORDER_POSITION
    ) => {
      setOrderedIDs((currentOrderedIDs) =>
        buildOrderedIDs(
          currentOrderedIDs,
          datasetID,
          targetDatasetID,
          orderPosition
        )
      );
    },
    []
  );

  // Saves order.
  const onSaveReorder = useCallback(() => {
    setMode(REORDER_MODE.INACTIVE);
    // TODO(cc) set up query to save order.
  }, []);

  // Starts reorder mode.
  const onStartReorder = useCallback(() => {
    setMode(REORDER_MODE.ACTIVE);
  }, []);

  // Sets initial order.
  useEffect(() => {
    if (!initialOrderedIDs) return;
    setOrderedIDs(initialOrderedIDs);
  }, [initialOrderedIDs]);

  return {
    isReorderUX,
    mode,
    reorderAction: {
      onCancelReorder,
      onReorder,
      onSaveReorder,
      onStartReorder,
    },
    orderedIDs,
  };
}

/**
 * Returns the updated order of datasets.
 * @param orderedIDs - Current dataset IDs, ordered.
 * @param datasetID - Dataset ID to reorder.
 * @param targetDatasetID - Target dataset ID to reorder to.
 * @param orderPosition - Indicates whether to insert before or after target dataset.
 * @returns order of datasets.
 */
function buildOrderedIDs(
  orderedIDs: string[],
  datasetID: string,
  targetDatasetID: string,
  orderPosition: ORDER_POSITION
): string[] {
  if (datasetID === targetDatasetID) return orderedIDs;
  const index = orderedIDs.indexOf(datasetID);
  const targetIndex = orderedIDs.indexOf(targetDatasetID);
  const nextOrder = [];
  for (let i = 0; i < orderedIDs.length; i++) {
    // Skip the dataset to reorder.
    if (i === index) continue;
    // Insert the dataset at the target index.
    if (i === targetIndex) {
      const newOrderedIDs = [datasetID, orderedIDs[i]];
      if (orderPosition === ORDER_POSITION.AFTER) {
        // Reverse order; the dataset should be inserted after the target dataset.
        newOrderedIDs.reverse();
      }
      nextOrder.push(...newOrderedIDs);
      continue;
    }
    // Add the rest of the datasets.
    nextOrder.push(orderedIDs[i]);
  }
  return nextOrder;
}
