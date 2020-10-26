import { get as getLocalStorage } from "src/common/localStorage/get";
import { BOOLEAN, set as setLocalStorage } from "src/common/localStorage/set";
import { isSSR } from "../utils/isSSR";

const FEATURE_FLAG_PREFIX = "cxg-ff-";

// Checks both URL and localStorage
export function get(key: string): string | null {
  if (isSSR()) return null;

  const params = new URLSearchParams(window.location.search);

  const URLValue = params.get(key);

  if (URLValue !== null) {
    const URLValueAsBooleanString: BOOLEAN =
      URLValue === "true" ? BOOLEAN.TRUE : BOOLEAN.FALSE;

    setLocalStorage(FEATURE_FLAG_PREFIX + key, URLValueAsBooleanString);
  }

  return getLocalStorage(FEATURE_FLAG_PREFIX + key);
}
