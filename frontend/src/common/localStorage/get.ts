import { isSSR } from "../utils/isSSR";

export function get(key: string): string | null {
  if (isSSR()) return null;

  return window?.localStorage?.getItem(key) || null;
}
