export const KEYS = {
  WORK_IN_PROGRESS_WARN: "cxg.WORK_IN_PROGRESS_WARN",
};

export function storageGet(
  key: string,
  defaultValue: string | null = null,
): string | null {
  try {
    const val = window.localStorage.getItem(key);
    if (val === null) return defaultValue;
    return val;
  } catch (e) {
    return defaultValue;
  }
}

/**
 * Removes local storage item.
 * @param key - local storage key name
 */
export function storageRemove(key: string): void {
  try {
    window.localStorage.removeItem(key);
  } catch {
    // continue
  }
}

export function storageSet(key: string, value: string): void {
  try {
    window.localStorage.setItem(key, value);
  } catch {
    // continue
  }
}
