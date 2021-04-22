import { get as getLocalStorage } from "src/common/localStorage/get";
import { BOOLEAN, set as setLocalStorage } from "src/common/localStorage/set";
import { isSSR } from "../utils/isSSR";
import { FEATURES } from "./features";

const FEATURE_FLAG_PREFIX = "cxg-ff-";

const allowedKeys = Object.values(FEATURES) as string[];

export function get(key: string): string | null {
  if (isSSR()) return null;

  return getLocalStorage(FEATURE_FLAG_PREFIX + key);
}

export function checkFeatureFlags(): void {
  if (isSSR()) return;

  const search = window.location.search;

  if (!search) return;

  const params = new URLSearchParams(search);

  params.forEach((value, key) => {
    if (!allowedKeys.includes(key)) return;

    setFeatureFlag(key, value);
  });
}

function setFeatureFlag(key: string, value: string) {
  const URLValueAsBooleanString =
    value === "true" ? BOOLEAN.TRUE : BOOLEAN.FALSE;

  setLocalStorage(FEATURE_FLAG_PREFIX + key, URLValueAsBooleanString);
}
