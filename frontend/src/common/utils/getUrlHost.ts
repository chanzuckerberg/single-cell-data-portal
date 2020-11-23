export function getUrlHost(url: string): string | null {
  try {
    const result = new URL(url);
    return result.host;
  } catch {
    return null;
  }
}
