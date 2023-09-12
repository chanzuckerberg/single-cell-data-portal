// Core dependencies
import { Intent } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { useRouter } from "next/router";
import { useEffect } from "react";
import { removeParams } from "src/common/utils/removeParams";
// App dependencies
import Toast from "src/views/Collection/components/Toast";

/**
 * Pop toast if user has been redirected from tombstoned collection.
 */
export function useExplainTombstoned(): void {
  const router = useRouter();
  const { tombstoned_collection_id } = router.query;
  const { isReady } = router;

  /* Check for tombstoned_collection_id query string param and pop toast if present. */
  useEffect(() => {
    if (isReady && tombstoned_collection_id) {
      Toast.show({
        icon: IconNames.ISSUE,
        intent: Intent.PRIMARY,
        message:
          "This collection was withdrawn. You’ve been redirected to Collections.",
      });
      removeParams("tombstoned_collection_id");
    }
  }, [tombstoned_collection_id, isReady]);
}
