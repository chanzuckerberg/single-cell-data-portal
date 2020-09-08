export function isSSR(): boolean {
  return typeof window === "undefined";
}
