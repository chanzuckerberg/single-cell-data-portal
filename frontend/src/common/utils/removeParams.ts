export function removeParams(params: Array<string> | string): void {
  if (typeof params === "string") params = [params];
  if (params.length < 1) return;
  const urlParams = new URLSearchParams(window.location.search);
  const url = window.location.href;
  const afterSlashBeforeParam = url
    .substring(url.indexOf("/") + 1)
    .split("?")[0];
  params.forEach((param) => urlParams.delete(param));
  let newURL = afterSlashBeforeParam;
  if (urlParams.toString().length > 0) {
    newURL += "?" + urlParams.toString();
  }
  window.history.replaceState(null, " ", "/" + newURL);
}
