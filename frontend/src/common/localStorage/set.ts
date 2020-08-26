export enum BOOLEAN {
  TRUE = "true",
  FALSE = "false",
}

export function set(key: string, value: BOOLEAN) {
  if (typeof window === "undefined") return;

  window?.localStorage?.setItem(key, value);
}
