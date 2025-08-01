export enum ROUTES {
  HOMEPAGE = "/",
  CENSUS_DIRECTORY = "/census-models",
  COLLECTION = "/collections/:id",
  COLLECTIONS = "/collections",
  DATASETS = "/datasets",
  DMCA = "/dmca/",
  TOS = "/tos/",
  PRIVACY = "/privacy/",
  PREVIEW_POLICIES = "/previewpolicies/",
  WHERE_IS_MY_GENE = "/gene-expression",
  DOCS = "/docs",
  PUBLISHED_DATA_DOCS = "/docs/02__Find%20Published%20Data",
  WMG_DOCS = "/docs/04__Analyze%20Public%20Data/4_2__Gene%20Expression%20Documentation/4_2_1__Get%20Started",
  FMG_DOCS = "/docs/04__Analyze%20Public%20Data/4_2__Gene%20Expression%20Documentation/4_2_5__Find%20Marker%20Genes",
  WMG_DOCS_ORDERING = "/docs/04__Analyze%20Public%20Data/4_2__Gene%20Expression%20Documentation/4_2_2__Cell%20Type%20and%20Gene%20Ordering",
  WMG_DOCS_DATA_PROCESSING = "/docs/04__Analyze%20Public%20Data/4_2__Gene%20Expression%20Documentation/4_2_3__Gene%20Expression%20Data%20Processing",
  SITEMAP = "/sitemap",
  DE = "/differential-expression",
  CELL_GUIDE = "/cellguide",
  CELL_GUIDE_CELL_TYPE = "/cellguide/:cellTypeId",
  CELL_GUIDE_TISSUE = "/cellguide/tissues/:tissueId",
  CELL_GUIDE_TISSUE_SPECIFIC_CELL_TYPE = "/cellguide/tissues/:tissueId/cell-types/:cellTypeId",
  DEPLOYED_VERSION = "/api/deployed_version",
}
