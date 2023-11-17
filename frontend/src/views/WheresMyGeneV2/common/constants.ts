export const GENE_SEARCH_BAR_HEIGHT_PX = 32;
export const HEATMAP_CONTAINER_ID = "heatmap-container-id";
export const FMG_GENE_STRENGTH_THRESHOLD = 0.5;

/**
 * (thuang): The `id` options here must match what the BE expects
 */
export const COMPARE_OPTIONS: {
  id: string | undefined;
  name: string;
}[] = [
  { id: undefined, name: "None" },
  { id: "disease", name: "Disease" },
  { id: "self_reported_ethnicity", name: "Ethnicity" },
  { id: "publication", name: "Publication" },
  { id: "sex", name: "Sex" },
];

export type CompareId = "disease" | "ethnicity" | "sex" | "publication";

export const getCompareOptionNameById = (id?: CompareId) => {
  const option = COMPARE_OPTIONS.find((option) => option.id === id);

  return option?.name || "None";
};

// This is used as the default/initial height of the x-axis
export const X_AXIS_CHART_HEIGHT_PX = 80;
export const X_AXIS_CHART_HEIGHT_PX_SVG = 30;

// This is the height of the container for the gene delete and hover icons
// Increasing this value adds more space between the gene label and icons
export const X_AXIS_HOVER_CONTAINER_HEIGHT_PX = 40;

export const MARGIN_BETWEEN_HEATMAPS = 4;

// Below constants are for left sidebar tooltip text so that tests can import these strings from this file instead
export const SELECT_TISSUE_GENE_TEXT =
  "Please select at least one tissue and gene to use this option.";
export const GROUP_BY_TOOLTIP_TEXT =
  "View expression for each cell type by the dimension selected.";
export const COLOR_SCALE_TOOLTIP_TEXT =
  "Expression is scaled to the range [0,1]. Scaling is done by assigning the minimum value in the current view to 0 and the max is assigned to 1.";
export const SORT_CELL_TYPES_TOOLTIP_TEXT =
  "Sort cell types by Cell Ontology or Hierarchical ordering. Cell ontology ordering groups cell types together based on their ontological relationships. Hierarchical ordering groups cell types with similar expression patterns together based on the genes selected.";
export const SORT_GENES_TOOLTIP_TEXT =
  "Sort genes As Entered or using Hierarchical ordering. Genes are displayed in the order they are added to the dot plot using As Entered ordering. Hierarchical ordering groups genes with similar expression patterns together.";

export const HOVER_START_TIME_MS = 2 * 1000;
