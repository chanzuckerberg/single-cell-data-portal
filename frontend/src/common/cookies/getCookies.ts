export function getCookie(key: string) {
  const name = key + "=";
  const decodedCookie = decodeURIComponent(document.cookie);
  const entries = decodedCookie.split(";");

  for (let i = 0; i < entries.length; i++) {
    let entry = entries[i];
    while (entry.charAt(0) === " ") {
      entry = entry.substring(1);
    }
    if (entry.indexOf(name) === 0) {
      return entry.substring(name.length, entry.length);
    }
  }

  return "";
}
