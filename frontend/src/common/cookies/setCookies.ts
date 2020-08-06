export enum BOOLEAN {
  TRUE = "true",
  FALSE = "false",
}

export function setCookie(name: string, value: BOOLEAN) {
  if (typeof document === "undefined") return false;

  const date = new Date();
  date.setMonth(date.getMonth() + 12);
  const nameString = `${name}=${value}; `;
  const expiryString = `expires=${date.toUTCString()}; `;
  const pathString = `path=/; `;
  const domainString =
    process.env.ENV !== "dev" ? `domain=${window.location.host}; ` : "";

  document.cookie = nameString + expiryString + pathString + domainString;
}
