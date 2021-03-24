export enum BOOLEAN {
  TRUE = "yes",
  FALSE = "no",
}

export function set(key: string, value: BOOLEAN) {
  if (typeof window === "undefined") return;

  window?.localStorage?.setItem(key, value);
}
