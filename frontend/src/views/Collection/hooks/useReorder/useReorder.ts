import { useCallback, useState } from "react";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { FEATURES } from "src/common/featureFlags/features";
import { useOrderDatasets } from "src/common/queries/collections";
import { Collection } from "src/common/entities";

export type OnCancelReorderFn = () => void;

export type OnReorderFn = (datasetIndex: number, targetIndex: number) => void;

export type OnSaveReorderFn = () => void;

export type OnStartReorderFn = (datasetIds: string[]) => void;

export interface ReorderAction {
  onCancelReorder: OnCancelReorderFn;
  onReorder: OnReorderFn;
  onSaveReorder: OnSaveReorderFn;
  onStartReorder: OnStartReorderFn;
}

export enum REORDER_MODE {
  ACTIVE = "ACTIVE",
  INACTIVE = "INACTIVE",
}

export interface UseReorder {
  isReorder: boolean;
  isReorderUX: boolean;
  orderedIds?: string[];
  reorderAction: ReorderAction;
}

/**
 * Reorder functionality for collection datasets.
 * The reorder mode can be either "inactive" or "active" and is used to enable or disable the datasets reorder feature
 * in the collection view.
 * @param collectionId - Collection ID to reorder datasets for.
 * @returns reorder mode.
 */
export function useReorder(collectionId: Collection["id"]): UseReorder {
  const isReorderUX = useFeatureFlag(FEATURES.REORDER); // Reorder datasets UX feature flag (reordering is currently only available with the feature flag).
  const [mode, setMode] = useState<REORDER_MODE>(REORDER_MODE.INACTIVE);
  const [orderedIds, setOrderedIds] = useState<string[]>();
  const orderDatasetsMutation = useOrderDatasets(collectionId);

  // Cancels reorder mode.
  const onCancelReorder = useCallback(() => {
    setMode(REORDER_MODE.INACTIVE);
    setOrderedIds(undefined);
  }, []);

  // Updates order.
  const onReorder = useCallback((datasetIndex: number, targetIndex: number) => {
    setOrderedIds((currentOrderedIds) =>
      buildOrderedIds(currentOrderedIds, datasetIndex, targetIndex)
    );
  }, []);

  // Saves order.
  const onSaveReorder = async (): Promise<void> => {
    setMode(REORDER_MODE.INACTIVE);

    // Send order to BE.
    const payload = JSON.stringify({
      datasets: orderedIds,
    });
    await orderDatasetsMutation.mutateAsync({
      collectionId,
      payload,
    });
  };

  // Starts reorder mode.
  const onStartReorder = useCallback((datasetIds: string[]) => {
    setMode(REORDER_MODE.ACTIVE);
    setOrderedIds(datasetIds);
  }, []);

  return {
    isReorder: mode === REORDER_MODE.ACTIVE,
    isReorderUX,
    reorderAction: {
      onCancelReorder,
      onReorder,
      onSaveReorder,
      onStartReorder,
    },
    orderedIds,
  };
}

/**
 * Returns the updated order of datasets.
 * @param orderedIds - Current dataset IDs, ordered.
 * @param datasetIndex - Index of dataset to reorder.
 * @param targetIndex - Index of target position to reorder dataset to.
 * @returns order of datasets.
 */
function buildOrderedIds(
  orderedIds: string[] | undefined,
  datasetIndex: number,
  targetIndex: number
): string[] | undefined {
  if (!orderedIds) return;
  // Reordering to the same position.
  if (datasetIndex === targetIndex) return orderedIds;
  const nextOrder = [...orderedIds];
  // Grab the dataset ID.
  const datasetId = nextOrder[datasetIndex];
  // Remove the dataset to reorder.
  nextOrder.splice(datasetIndex, 1);
  // Insert the dataset at the target index.
  nextOrder.splice(targetIndex, 0, datasetId);
  return nextOrder;
}
