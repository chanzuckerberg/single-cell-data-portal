export function get(key: string): string | null {
  return window?.localStorage?.getItem(key) || null;
}
