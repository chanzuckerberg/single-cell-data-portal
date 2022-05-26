import { useEffect, useState } from "react";
import { getObject } from "src/common/sessionStorage/get";
import { KEYS, setObject } from "src/common/sessionStorage/set";

/**
 * Save value to session storage on change. Initial value is read from session storage and if it's not yet set, it's
 * set to the given default value.
 * @param key - Name of key to use for session storage.
 * @param defaultValue - Value to use if session storage with the given key is not yet set.
 * @returns Value saved in session storage (or default value), and mutator for setting hook's internal value.
 */
export function useSessionStorage<T>(
  key: KEYS,
  defaultValue: T
): [T, (nextValue: T) => void] {
  // Value to save to session storage. Initial state read from session storage, otherwise set to default value.
  const [value, setValue] = useState<T>(() => {
    const savedValue = getObject<T>(key);
    // Explicit null check as T can possibly be boolean.
    return savedValue === null ? defaultValue : savedValue;
  });

  // Save value to sessions storage on change.
  useEffect(() => {
    setObject(key, value);
  }, [key, value]);

  return [value, setValue];
}
