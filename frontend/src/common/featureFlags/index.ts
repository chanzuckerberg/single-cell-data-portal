import { get as getLocalStorage } from "src/common/localStorage/get";
import { BOOLEAN, set as setLocalStorage } from "src/common/localStorage/set";
import { isSSR } from "../utils/isSSR";
import { removeParams } from "../utils/removeParams";
import { FEATURES } from "./features";

const FEATURE_FLAG_PREFIX = "cxg-ff-";

// Checks both URL and localStorage
export function get(key: string): string | null {
  if (isSSR()) return null;

  const newlyAddedParams: Array<string> = [];
  const params = new URLSearchParams(window.location.search);
  if (window.location.search.indexOf("?") !== -1) {
    Object.values(FEATURES).forEach((feature) => {
      if (params.has(feature)) {
        const URLValueAsBooleanString =
          params.get(feature) === "true" ? BOOLEAN.TRUE : BOOLEAN.FALSE;
        setLocalStorage(FEATURE_FLAG_PREFIX + feature, URLValueAsBooleanString);
        newlyAddedParams.push(feature);
      }
    });
  }

  removeParams(newlyAddedParams);

  return getLocalStorage(FEATURE_FLAG_PREFIX + key);
}
