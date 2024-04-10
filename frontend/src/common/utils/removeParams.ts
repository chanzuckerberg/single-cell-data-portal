import { NextRouter } from "next/router";

export function removeParams({
  params,
  router,
}: {
  params: Array<string> | string;
  router: NextRouter;
}): void {
  if (!params || !params.length) return;

  if (typeof params === "string") params = [params];

  const { search, href } = window.location;

  const urlParams = new URLSearchParams(search);
  const pathname = new URL(href).pathname;

  params.forEach((param) => urlParams.delete(param));

  const query = urlParams.toString();
  const newURL = pathname + (query ? "?" + query : "");

  router.replace(newURL);
}
