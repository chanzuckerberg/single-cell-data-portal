import { QueryGroup } from "src/views/DifferentialExpression/common/store/reducer";
import { LATEST_SHARE_LINK_VERSION } from "src/views/WheresMyGeneV2/components/GeneSearchBar/components/ShareButton/utils";

const HUMAN_ORGANISM_ID = "NCBITaxon:9606";

export const generateAndCopyShareUrl = ({
  queryGroup,
  organism,
  genes,
  cellTypes = [],
}: {
  queryGroup: QueryGroup;
  organism: string;
  genes: string[];
  cellTypes: string[];
}) => {
  const url = new URL(window.location.href);
  url.pathname = "gene-expression";
  // human is empty default
  if (organism !== HUMAN_ORGANISM_ID) {
    url.searchParams.set("organism", organism);
  }
  const filters: Partial<Record<string, string[]>> = {};
  if (queryGroup.diseases.length > 0) filters.diseases = queryGroup.diseases;
  if (queryGroup.ethnicities.length > 0)
    filters.ethnicities = queryGroup.ethnicities;
  if (queryGroup.publicationCitations.length > 0)
    filters.publications = queryGroup.publicationCitations;
  if (queryGroup.sexes.length > 0) filters.sexes = queryGroup.sexes;
  if (queryGroup.tissues.length > 0) filters.tissues = queryGroup.tissues;

  Object.entries(filters).forEach(([key, value]) => {
    value && url.searchParams.set(key, value.join(","));
  });
  url.searchParams.set("genes", genes.join(","));

  if (cellTypes.length > 0) {
    url.searchParams.set("cellTypes", cellTypes.join(","));
  }

  url.searchParams.set("ver", LATEST_SHARE_LINK_VERSION);

  return String(url);
};
