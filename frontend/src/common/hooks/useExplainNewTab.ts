// Core dependencies
import { Intent } from "@blueprintjs/core";
import { useRouter } from "next/router";
import { useEffect } from "react";
// App dependencies
import Toast from "src/views/Collection/components/Toast";

/**
 * Pop toast with given message if user has come from Explorer with work in progress.
 * @param message - Text to display in toast.
 */
export function useExplainNewTab(message: string): void {
  const router = useRouter();
  const { explainNewTab } = router.query;
  const { isReady } = router;

  /* Check for explainNewTab query string param and pop toast if present. */
  useEffect(() => {
    // Query param with no value (e.g. explainNewTab in http://url.com/?explainNewTab) resolves to empty string.
    if (isReady && explainNewTab === "") {
      Toast.show({
        intent: Intent.PRIMARY,
        message,
      });
    }
  }, [explainNewTab, isReady, message]);
}
