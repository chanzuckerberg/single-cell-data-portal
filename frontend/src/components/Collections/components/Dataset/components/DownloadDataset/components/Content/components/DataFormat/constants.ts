import { DATASET_ASSET_FORMAT } from "src/common/entities";
export const POSSIBLE_DOWNLOAD_FORMATS = [
  {
    format: DATASET_ASSET_FORMAT.H5AD,
    label: ".h5ad (AnnData v0.10)",
    type: "RNA",
    description:
      "Hierarchical data format used to store annotated data matrices, typically for single-cell omics data, including gene expression, metadata, and embeddings.",
  },
  {
    format: DATASET_ASSET_FORMAT.ATAC_INDEX,
    label: ".tsv (Fragments w/ index)",
    type: "DNA ACCESSIBILITY",
    description:
      "Tab-separated values file listing sequencing ATAC-seq fragments with an accompanying index for fast genomic range queries.",
  },
];

export const getNotAvailableText = (filetype: string) =>
  `${filetype} file is unavailable for download at the moment`;
