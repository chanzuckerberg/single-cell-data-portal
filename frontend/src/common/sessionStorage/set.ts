/**
 * Key names of valid session storage keys.
 */
export enum KEYS {
  "FILTER_DATASETS" = "cxg-filter-datasets",
  "FILTER_DATASETS_SELECTED_UI" = "cxg-filter-datasets-selected-ui",
  "FILTER_COLLECTIONS" = "cxg-filter-collections",
  "FILTER_COLLECTIONS_SELECTED_UI" = "cxg-filter-collections-selected-ui",
  "SIDE_BAR_DATASETS" = "cxg-sidebar-datasets",
  "SIDE_BAR_COLLECTIONS" = "cxg-sidebar-collections",
}

/**
 * Add string value to session storage.
 * @param key - Key name of value to add to session storage.
 * @param value - String value to add to session storage.
 */
export function set(key: KEYS, value: string): void {
  if (typeof window === "undefined") return;
  window?.sessionStorage?.setItem(key, value);
}

/**
 * Add object to session storage.
 * @param key - Key name of value to add to session storage.
 * @param value - Object value to add to session storage.
 */
export function setObject<T>(key: KEYS, value: T): void {
  if (typeof window === "undefined") return;
  window?.sessionStorage?.setItem(key, JSON.stringify(value));
}
