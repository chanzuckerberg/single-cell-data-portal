// App dependencies
import { KEYS } from "src/common/sessionStorage/set";
import { isSSR } from "../utils/isSSR";

/**
 * Return the session storage string value with the given key.
 * @param key - Key name of session storage value to return.
 * @returns Stored (string) value for the given session storage key.
 */
export function get(key: KEYS): string | null {
  if (isSSR()) return null;

  return window?.sessionStorage?.getItem(key) ?? null;
}

/**
 * Return the session storage (object) value with the given key.
 * @param key - Key name of session storage value to return.
 * @returns Stored (object) value for the given session storage key.
 */
export function getObject<T>(key: KEYS): T | null {
  if (isSSR()) return null;

  const value = window?.sessionStorage?.getItem(key);
  if (!value) {
    return null;
  }

  // Attempt to convert storage value to object.
  try {
    return JSON.parse(value);
  } catch (e) {
    // Return null if object can't be parsed.
    return null;
  }
}
