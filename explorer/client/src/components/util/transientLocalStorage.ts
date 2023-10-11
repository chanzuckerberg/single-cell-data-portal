/* App dependencies */
import { storageGet, storageRemove, storageSet } from "./localStorage";

/**
 * Gets local storage that has been set with an expiry value. Single use as local storage is cleared on get.
 * @param key - local storage key name
 * @returns True if expiry time has not expired.
 */
export function storageGetTransient(key: string): boolean {
  const storeValue = storageGet(key);
  storageRemove(key);

  try {
    if (storeValue) {
      const expiryTime = JSON.parse(storeValue);
      return new Date().getTime() < expiryTime;
    }
    return false;
  } catch (e) {
    return false;
  }
}

/**
 * Sets local storage with an expiry value.
 * @param key - local storage key name
 * @param timeout - time in ms
 */
export function storageSetTransient(key: string, timeout: number): void {
  const now = new Date();
  const expiryTime = now.getTime() + timeout;
  storageSet(key, JSON.stringify(expiryTime));
}
