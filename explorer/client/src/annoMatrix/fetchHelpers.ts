export { doBinaryRequest, doFetch } from "../util/actionHelpers";

/* double URI encode - needed for query-param filters */
export function _dubEncURIComp(s: string | number | boolean): string {
  return encodeURIComponent(encodeURIComponent(s));
}
