export function get(key: string): string | null {
  if (typeof window === "undefined") return null;

  return window?.localStorage?.getItem(key) || null;
}
