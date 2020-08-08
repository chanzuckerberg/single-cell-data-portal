export enum BOOLEAN {
  TRUE = "true",
  FALSE = "false",
}

export function set(key: string, value: BOOLEAN) {
  window?.localStorage?.setItem(key, value);
}
