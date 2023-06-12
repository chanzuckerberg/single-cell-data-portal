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

// This is used as the default/initial height of the x-axis
export const X_AXIS_CHART_HEIGHT_PX = 80;

// This is the height of the container for the gene delete and hover icons
// Increasing this value adds more space between the gene label and icons
export const X_AXIS_HOVER_CONTAINER_HEIGHT_PX = 40;

export const MARGIN_BETWEEN_HEATMAPS = 8;
