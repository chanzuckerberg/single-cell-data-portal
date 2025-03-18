import { DATASET_ASSET_FORMAT } from "src/common/entities";
export const possibleDownloadFormats = [
  {
    format: DATASET_ASSET_FORMAT.H5AD,
    label: ".h5ad (AnnData v0.10)",
    type: "RNA",
    description:
      "Lorem ipsum odor amet, consectetuer adipiscing elit. Augue nulla etiam dolor orci praesent.",
  },
  {
    format: DATASET_ASSET_FORMAT.ATAC_INDEX,
    label: ".tsv (Fragments w/ index)",
    type: "DNA ACCESSIBILITY",
    description:
      "Turpis habitasse potenti dapibus tincidunt fames metus vulputate feugiat.",
  },
];
