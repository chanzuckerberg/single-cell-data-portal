import { useCallback, useState } from "react";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { FEATURES } from "src/common/featureFlags/features";
import { useOrderDatasets } from "src/common/queries/collections";
import { Collection } from "src/common/entities";

export type OnCancelReorderFn = () => void;

export type OnReorderFn = (
  datasetID: string,
  targetDatasetID: string,
  position: ORDER_POSITION
) => void;

export type OnSaveReorderFn = () => void;

export type OnStartReorderFn = (datasetIDs: string[]) => void;

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
  orderedIDs?: string[];
  reorderAction: ReorderAction;
}

/**
 * Reorder functionality for collection datasets.
 * The reorder mode can be either "inactive" or "active" and is used to enable or disable the datasets reorder feature
 * in the collection view.
 * @param collectionId - ID of collection to reorder datasets for.
 * @returns reorder mode.
 */
export function useReorderMode(collectionId: Collection["id"]): UseReorderMode {
  const isReorderUX = useFeatureFlag(FEATURES.REORDER); // Reorder datasets UX feature flag (reordering is currently only available with the feature flag).
  const [mode, setMode] = useState<REORDER_MODE>(REORDER_MODE.INACTIVE);
  const [orderedIDs, setOrderedIDs] = useState<string[]>();
  const orderDatasetsMutation = useOrderDatasets(collectionId);

  // Cancels reorder mode.
  const onCancelReorder = useCallback(() => {
    setMode(REORDER_MODE.INACTIVE);
    setOrderedIDs(undefined);
  }, []);

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
  const onSaveReorder = async (): Promise<void> => {
    setMode(REORDER_MODE.INACTIVE);

    // Send order to BE.
    const payload = JSON.stringify({
      datasets: orderedIDs,
    });
    await orderDatasetsMutation.mutateAsync({
      collectionId,
      payload,
    });
  };

  // Starts reorder mode.
  const onStartReorder = useCallback((datasetIDs: string[]) => {
    setMode(REORDER_MODE.ACTIVE);
    setOrderedIDs(datasetIDs);
  }, []);

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
  orderedIDs: string[] | undefined,
  datasetID: string,
  targetDatasetID: string,
  orderPosition: ORDER_POSITION
): string[] | undefined {
  if (!orderedIDs) return;
  // Reordering to the same position.
  if (datasetID === targetDatasetID) return orderedIDs;
  const nextOrder = [...orderedIDs];
  // Remove the dataset to reorder.
  const index = nextOrder.indexOf(datasetID);
  nextOrder.splice(index, 1);
  // Remove the target dataset.
  const targetIndex = nextOrder.indexOf(targetDatasetID);
  nextOrder.splice(targetIndex, 1);
  // Insert the dataset and target dataset at the target index, in the correct order.
  const datasetIDs = [datasetID, targetDatasetID];
  if (orderPosition === ORDER_POSITION.AFTER) {
    // Reverse order; the dataset should be inserted after the target dataset.
    datasetIDs.reverse();
  }
  nextOrder.splice(targetIndex, 0, ...datasetIDs);
  return nextOrder;
}
