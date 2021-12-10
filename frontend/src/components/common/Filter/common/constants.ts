/* Case insensitive sorter */
export const COLLATOR_CASE_INSENSITIVE = new Intl.Collator("en", {
  numeric: true,
  sensitivity: "base",
});
