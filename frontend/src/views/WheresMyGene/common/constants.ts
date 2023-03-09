export const HEATMAP_CONTAINER_ID = "heatmap-container-id";
export const FMG_GENE_STRENGTH_THRESHOLD = 0.5;

/**
 * (thuang): The `id` options here must match what the BE expects
 */
export const COMPARE_OPTIONS = [
  { id: undefined, name: "None" },
  { id: "disease", name: "Disease" },
  { id: "self_reported_ethnicity", name: "Ethnicity" },
  { id: "sex", name: "Sex" },
];

export type CompareId = "disease" | "ethnicity" | "sex";

export const getCompareOptionNameById = (id?: CompareId) => {
  const option = COMPARE_OPTIONS.find((option) => option.id === id);

  return option?.name || "None";
};
