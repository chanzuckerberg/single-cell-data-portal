export enum ROUTES {
  HOMEPAGE = "/",
  MY_COLLECTIONS = "/my-collections",
  COLLECTION = "/collections/:id",
  COLLECTIONS = "/collections",
  DATASETS = "/datasets",
  TOS = "/tos/",
  PRIVACY = "/privacy/",
  PREVIEW_POLICIES = "/previewpolicies/",
  WHERE_IS_MY_GENE = "/scExpression",
}

export enum EXTERNAL_LINKS {
  // (thuang): TEMP -- revert back to the old links once Google Analytics is
  // removed from Gitbook
  // DOCS_DATA_PORTAL = "https://docs.cellxgene.cziscience.com",
  // DOCS_ROADMAP = "https://docs.cellxgene.cziscience.com/roadmap",
  // DOCS_TUTORIAL = "https://docs.cellxgene.cziscience.com/portal/data-portal",
  DOCS_DATA_PORTAL = "https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/README.md",
  DOCS_ROADMAP = "https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/roadmap.md",
  DOCS_TUTORIAL = "https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/portal/data-portal.md",
  WMG_DOC = "https://github.com/chanzuckerberg/cellxgene-documentation/blob/pablo-gar/wheres-my-gene/scExpression/scExpression-documentation.md",
}
